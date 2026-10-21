import hashlib
import json
from pathlib import Path
import subprocess
import sys

from fullrun_common import ROOT, OUT, BASELINE_COMMIT, configure, verify_source, execution_commit
configure()


def main():
    import numpy as np
    from scipy.stats import kstest,spearmanr
    import matplotlib.pyplot as plt
    from colossus import settings
    settings.BASE_DIR=str(OUT/'analysis_cache')
    sys.path.insert(0,str(ROOT/'Insert_Halo_v1.00'))
    from IH_ABC_strong_lensing import IH_ABC
    rows=[];batches=[]
    for i in range(3):
        r=json.loads((OUT/f'abc_shard_{i}'/'abc_full.json').read_text())
        batches.append(r);rows.extend([dict(t,shard=i) for t in r['trials']])
    assert len(rows)==200 and all(r['status']=='completed' for r in batches)
    theta=np.array([[r['parameters']['sigma_0'],r['parameters']['vel_scale']] for r in rows])
    logs=np.log10(theta)
    scores=np.array([r['S_all_references'] for r in rows])
    flux=np.array([r['fluxes'] for r in rows])
    assert np.isfinite(theta).all() and np.isfinite(scores).all() and np.isfinite(flux).all()
    assert len(np.unique(theta,axis=0))==200
    assert ((theta[:,0]>=1e-3)&(theta[:,0]<=100)&(theta[:,1]>=1)&(theta[:,1]<=1000)).all()
    baseline=np.array(batches[0]['baseline_fluxes'])
    assert all(np.array_equal(baseline,np.array(b['baseline_fluxes'])) for b in batches)
    assert all(b['obs_fluxes']==batches[0]['obs_fluxes'] for b in batches)
    sorter=IH_ABC()
    sorter.drawn_parameter_values=np.array([r['parameters'] for r in rows],dtype=object)
    sorter.summary_statistics=scores[:,0]
    chosen,_,_=sorter.sort_sum_stats(acc_fraction=.2)
    assert len(chosen)==40
    by_reference={}
    for j,name in enumerate(['B','C','D','A']):
        sorter.summary_statistics=scores[:,j]
        idx,_,_=sorter.sort_sum_stats(acc_fraction=.2)
        by_reference[name]={'accepted_indices':idx.tolist(),'epsilon':float(scores[idx[-1],j]),
            'S_percentiles':np.percentile(scores[:,j],[0,5,20,50,80,95,100]).tolist(),
            'overlap_with_B':int(len(set(idx)&set(chosen)))}
    tests={}
    for j,(name,lo,hi) in enumerate([('sigma_0',-3,2),('vel_scale',0,3)]):
        ks=kstest((logs[chosen,j]-lo)/(hi-lo),'uniform')
        corr=spearmanr(logs[:,j],scores[:,0])
        counts,edges=np.histogram(logs[chosen,j],bins=np.linspace(lo,hi,21))
        tests[name]={'accepted_log10_histogram':counts.tolist(),'histogram_edges':edges.tolist(),
            'accepted_log10_quantiles':np.percentile(logs[chosen,j],[0,16,50,84,100]).tolist(),
            'KS_uniform_statistic':float(ks.statistic),'KS_uniform_pvalue':float(ks.pvalue),
            'spearman_score_statistic':float(corr.statistic),'spearman_score_pvalue':float(corr.pvalue)}
    maps=[]
    targets={1:.07164411564009043,10:.04709088674559325,100:.18557687608833362}
    for sigma in [1,10,100]:
        r=json.loads((OUT/f'sigma{sigma}_w100_seed20261021.json').read_text())
        assert r['status']=='completed'
        r['paper_kappa_max']=targets[sigma]
        r['kappa_max_ratio_to_paper']=r['kappa_max']/targets[sigma]
        maps.append(r)
    result={'main_commit':BASELINE_COMMIT, 'driver_commit':execution_commit(),
        'abc_completed':200,'accepted':40,'reference_image':'B','acceptance_fraction':.2,
        'baseline_S_by_B_C_D_A':batches[0]['baseline_S_all_references'],
        'observed_fluxes_B_C_D_A':batches[0]['obs_fluxes'],
        'baseline_fluxes_B_C_D_A':baseline.tolist(),'reference_sensitivity':by_reference,
        'parameter_statistics':tests,'unique_parameter_draws':len(np.unique(theta,axis=0)),
        'accepted_parameter_draws':theta[chosen].tolist(),
        'sum_trial_seconds':sum(r['seconds'] for r in rows),
        'trial_seconds_percentiles':np.percentile([r['seconds'] for r in rows],[0,50,90,95,100]).tolist(),
        'peak_process_rss_MiB':max(r['max_rss_MiB'] for r in rows),
        'adjustment_fraction_percentiles':np.percentile([r['mass_adjustment_fraction'] for r in rows],[0,50,90,95,100]).tolist(),
        'trials_with_adjustment_over_5pct':sum(r['mass_adjustment_fraction']>.05 for r in rows),
        'fractional_flux_change_mean':np.mean(flux/baseline-1,axis=0).tolist(),
        'fractional_flux_change_std':np.std(flux/baseline-1,axis=0).tolist(),
        'fractional_flux_change_min':np.min(flux/baseline-1,axis=0).tolist(),
        'fractional_flux_change_max':np.max(flux/baseline-1,axis=0).tolist(),
        'subhalo_count_range':[min(r['subhalos'] for r in rows),max(r['subhalos'] for r in rows)],
        'fieldhalo_count_range':[min(r['fieldhalos'] for r in rows),max(r['fieldhalos'] for r in rows)],
        'lens_component_range':[min(r['lens_components'] for r in rows),max(r['lens_components'] for r in rows)],
        'convergence_maps':maps}
    independent=[]
    for sigma,seed in [(1,20261022),(10,20261023),(100,20261024)]:
        p=OUT/f'sigma{sigma}_w100_seed{seed}.json'
        if p.exists():
            item=json.loads(p.read_text());assert item['status']=='completed'
            item['paper_kappa_max']=targets[sigma]
            independent.append(item)
    assert len(independent)==3
    result['independent_seed_convergence_maps']=independent
    imgstats=OUT/'paper_image_statistics.json'
    if imgstats.exists():
        result['paper_image_statistics']=json.loads(imgstats.read_text())
        for item,sigma in zip(maps,[1,10,100]):
            ref=result['paper_image_statistics']['convergence_maps'][str(sigma)]['kappa_percentiles_approx']
            item['paper_kappa_percentiles_approx']=ref
            item['median_kappa_ratio_to_paper_approx']=item['kappa_percentiles'][0]/ref[0]
    controls=OUT/'baseline_controls.json'
    if not controls.is_file():
        raise RuntimeError('Complete baseline-controls outputs are required for full-run analysis')
    if controls.exists():
        result['baseline_controls']=json.loads(controls.read_text())
        assert len(result['baseline_controls'])==3
        bs=[batches[0]['baseline_S_all_references'][0]]+[r['baseline_S_all_references'][0] for r in result['baseline_controls']]
        result['baseline_PSO_seed_range_B']=[min(bs),max(bs)]
        result['baseline_PSO_seed_std_B']=float(np.std(bs))
        result['ABC_S_std_B']=float(np.std(scores[:,0]))
    for name in ['profile_accuracy_full','wgd_fixed_case']:
        result[name]=json.loads((OUT/f'{name}.json').read_text())
    diag=result['wgd_fixed_case']
    result['fixed_case_fractional_magnification_changes_A_B_C_D']=(
        np.asarray(diag['magnification_substructure'])/np.asarray(diag['magnification_baseline'])-1).tolist()
    original=verify_source()
    result['unchanged_research_files']=len(original['research_files_sha256'])
    (OUT/'fullrun_summary.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(OUT/'abc_all200.npz',theta=theta,scores_all_references=scores,fluxes=flux,
                        accepted_B=chosen)

    plt.rcParams.update({'font.size':11})
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for j,(ax,name,lo,hi) in enumerate(zip(axes[:2],[r'$\log_{10}(\sigma_0/m_\chi)$',r'$\log_{10}w$'],[-3,0],[2,3])):
        edges=np.linspace(lo,hi,21)
        ax.hist(logs[chosen,j],bins=edges,color='tab:blue',alpha=.75,edgecolor='white')
        ax.hist(logs[:,j],bins=edges,weights=np.full(200,.2),histtype='step',color='black',lw=1.,label='All 200 draws x 0.2')
        ax.axhline(2,color='gray',ls='--',label='Uniform-prior expectation (40/20)')
        ax.set(xlabel=name,ylabel='Accepted realizations',title='40 accepted of 200 full trials')
        ax.set_ylim(0,int(np.histogram(logs[chosen,j],bins=edges)[0].max())+2)
        ax.legend(fontsize=8,loc='upper right')
    sc=axes[2].scatter(logs[:,0],logs[:,1],c=scores[:,0],s=15,cmap='viridis')
    axes[2].scatter(logs[chosen,0],logs[chosen,1],facecolors='none',edgecolors='red',s=50,label='Best 20%')
    axes[2].set(xlabel=r'$\log_{10}(\sigma_0/m_\chi)$',ylabel=r'$\log_{10}w$')
    axes[2].legend(fontsize=8);fig.colorbar(sc,ax=axes[2],label='S (reference B)')
    fig.savefig(OUT/'abc_full_results.png',dpi=160);plt.close(fig)

    arrays=[np.load(OUT/f'sigma{s}_w100_seed20261021_kappa.npz')['kappa'] for s in [1,10,100]]
    vmin=-2.93;vmax=np.log10(max(.18557687608833362,max(float(k.max()) for k in arrays)))
    fig,axes=plt.subplots(1,3,figsize=(14,4.7),layout='constrained')
    for ax,k,s in zip(axes,arrays,[1,10,100]):
        im=ax.imshow(np.log10(k),origin='lower',extent=[-3,3,-3,3],cmap='turbo',vmin=vmin,vmax=vmax)
        ax.set(xlabel='arcsec',ylabel='arcsec',title=rf'$\sigma_0/m_\chi={s}$; $\kappa_{{max}}={k.max():.4f}$')
    fig.colorbar(im,ax=axes,label=r'$\log_{10}\kappa_{sub}$')
    fig.savefig(OUT/'convergence_full_comparison.png',dpi=160);plt.close(fig)
    if len(independent)==3:
        fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
        for row,cases in enumerate([maps,independent]):
            for ax,item in zip(axes[row],cases):
                cfg=item['configuration'];s=cfg['sigma_0'];seed=cfg['seed']
                k=np.load(OUT/f'sigma{s:g}_w100_seed{seed}_kappa.npz')['kappa']
                im=ax.imshow(np.log10(k),origin='lower',extent=[-3,3,-3,3],cmap='turbo',vmin=vmin,vmax=vmax)
                ax.set(xlabel='arcsec',ylabel='arcsec',title=rf'$\sigma_0/m_\chi={s:g}$; $\kappa_{{max}}={k.max():.4f}$')
                ax.text(.03,.97,f'seed {seed}',va='top',transform=ax.transAxes,fontsize=9,color='white')
        fig.colorbar(im,ax=axes,label=r'$\log_{10}\kappa_{sub}$')
        fig.savefig(OUT/'convergence_all_six.png',dpi=160);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
