import json
from pathlib import Path
import pickle
import sys
import time

from fullrun_common import OUT, SRC, configure
configure()


def profile_accuracy():
    import numpy as np
    from scipy.integrate import quad
    from IH_alternative_DM_dist import gaussian_decomposition
    from IH_halo_conversion import HaloConversion
    from fullrun_maps import run
    run(100.,25.,20261021,False)
    run(10.,24.33,20261021,False)
    a=np.load(OUT/'sigma10_w24.33_seed20261021_halos.npz')
    gd=gaussian_decomposition(8,50)
    details={};t0=time.perf_counter()
    for kind in ['sub','field']:
        errors=[];calibration_errors=[];physical_ratios=[];projection_ratios=[];negative_counts=[];all_records=[]
        for i in range(len(a[kind+'_mass'])):
            p={k:float(a[f'{kind}_{k}'][i]) for k in ['rs','rc','rhos']}
            if kind=='sub':
                p.update(profile='TNFWC_SH_tau',rt=float(a[kind+'_rt'][i]));fn=HaloConversion.TNFWC_SH_tau_f
            else:
                p.update(profile='TNFWC_SH');fn=HaloConversion.TNFWC_SH_f
            amp,width,_=gd.decompose_direct_los_all(np.log10(p['rc']*1e-3),np.log10(p['rs']*1e7),fn,p)
            direct=float(gd.lens_projected_direct_all(fn,p,np.array([1e-3]))[0])
            reconstructed=float(np.sum(amp[amp>0]*np.exp(-1e-6/(2*width[amp>0]**2))))
            ca,cw,_=gd.decompose_direct_los_all(np.log10(p['rc']*1e-3),np.log10(p['rs']*1e9),fn,p)
            reconstructed_calibration=float(np.sum(ca[ca>0]*np.exp(-1e-6/(2*cw[ca>0]**2))))
            projected_radius=1e-3/p['rs']
            if kind=='sub':
                integrand=lambda t:float(fn(np.sqrt(projected_radius**2+t*t),p['rt']/p['rs'],p['rc']/p['rs'],p['rhos']))
            else:
                integrand=lambda t:float(fn(np.sqrt(projected_radius**2+t*t),p['rc']/p['rs'],p['rhos']))
            physical=2*p['rs']*quad(integrand,0,np.inf,epsabs=1e-4,epsrel=1e-7,limit=200)[0]
            err=abs(reconstructed/direct-1)
            errors.append(err);calibration_errors.append(abs(reconstructed_calibration/direct-1))
            physical_ratios.append(reconstructed/physical);projection_ratios.append(direct/physical)
            negative_counts.append(int(np.sum(amp<=0)))
            all_records.append({'mass':float(a[kind+'_mass'][i]),'t_over_tc':float(a[kind+'_t_over_tc'][i]),
                                'relative_error':err,'profile':p})
            if i and i%500==0:print(json.dumps({'profile_check':kind,'done':i}),flush=True)
        details[kind]={'count':len(errors),'error_percentiles':np.percentile(errors,[0,50,90,95,99,100]).tolist(),
            'calibration_GDmax_1e9_error_percentiles':np.percentile(calibration_errors,[0,50,90,95,99,100]).tolist(),
            'production_GDmax_1e7_error_percentiles':np.percentile(errors,[0,50,90,95,99,100]).tolist(),
            'fraction_error_above_20pct':float(np.mean(np.asarray(errors)>.2)),
            'physical_normalization_ratio_percentiles':np.percentile(physical_ratios,[0,50,95,100]).tolist(),
            'code_projection_over_full_physical_percentiles':np.percentile(projection_ratios,[0,50,95,100]).tolist(),
            'maximum_discarded_gaussians':max(negative_counts),
            'worst':max(all_records,key=lambda r:r['relative_error'])}
    details['seconds']=time.perf_counter()-t0
    (OUT/'profile_accuracy_full.json').write_text(json.dumps(details,indent=2))
    print(json.dumps({'profile_accuracy_completed':details}),flush=True)


def wgd_fixed_case():
    import numpy as np
    from IH_main import Insert_Halo
    from IH_data_analysis import Insert_Halo_Data_Analysis
    from IH_ABC_strong_lensing import IH_ABC
    with (OUT/'fitted_wgd.pkl').open('rb') as f:s=pickle.load(f)
    s.setup_data_modelplot();obs,labels,_,_=s.get_flux_image_from_data(flux_threshold_log=.25)
    s.setup_lightmodels(add_lens_light=True,add_shapelets=True)
    s.setup_lensmodel(data_only=True);s.setup_image_model()
    kw=s.fitted_model_kwargs
    common=dict(kwargs_source=kw['kwargs_source'],kwargs_lens_light=kw['kwargs_lens_light'],kwargs_ps=kw['kwargs_ps'])
    base=s.image_model.image(kwargs_lens=kw['kwargs_lens'],**common)
    point_base=s.image_model.image(kwargs_lens=kw['kwargs_lens'],source_add=False,lens_light_add=False,**common)
    x=kw['kwargs_ps'][0]['ra_image'];y=kw['kwargs_ps'][0]['dec_image']
    xp,yp=s.image_model.Data.map_coord2pix(x,y)
    point_labels=[['B','C','D','A'][int(labels[int(round(py)),int(round(px))])-1]
                  if labels[int(round(py)),int(round(px))]>0 else 'outside flux mask'
                  for px,py in zip(xp,yp)]
    mu_base=s.usable_lensmodel.magnification(x,y,kw['kwargs_lens'])
    base_flux=np.array(s.get_flux_image_from_sim(base,labels))
    s.reset_image_model()
    cfg=dict(zlens=.23,zsource=.78,Mhost=13.5,msub_min=6,msub_max=11,cone_angle_arcsec=s.image_width_arcsec,
        center_x=s.x_center,center_y=s.y_center,DM_type='SIDM',sigma_0=100.,vel_scale=24.33,
        spatial='einasto',ngaussian=50,seed=20261021,fix_mass_range_sh=True,sigmalogc=.13,
        sh_dz=.01,sh_N_ma=500,sh_N_herm=20,sh_zmax=10.)
    print(json.dumps({'starting_wgd_fixed_case':cfg}),flush=True)
    t0=time.perf_counter();np.random.seed(20261021)
    ih=Insert_Halo();hc=ih.run_insert_halo(**cfg)
    da=Insert_Halo_Data_Analysis(hc)
    lm=da.prepare_lensmodel(mass_sheet_correction=False,include_subhalos=True,include_fieldhalos=True,
                            include_host=False,pre_model_return=True)
    s.setup_lensmodel(substructure_lensmodel_data=lm);s.setup_image_model()
    total_lens=kw['kwargs_lens']+lm['kwargs_substructure']
    full=s.image_model.image(kwargs_lens=total_lens,**common)
    point_full=s.image_model.image(kwargs_lens=total_lens,source_add=False,lens_light_add=False,**common)
    mu_full=s.usable_lensmodel.magnification(x,y,total_lens)
    full_flux=np.array(s.get_flux_image_from_sim(full,labels))
    abc=IH_ABC()
    score=[float(np.linalg.norm(np.array(abc.compute_flux_ratios(obs,i))-np.array(abc.compute_flux_ratios(full_flux,i)))) for i in range(4)]
    np.savez_compressed(OUT/'wgd_fixed_case_images.npz',baseline=base,substructure=full,
                        point_baseline=point_base,point_substructure=point_full,
                        observed=s.multi_band_list[0][0]['image_data'])
    r={'configuration':cfg,'seconds':time.perf_counter()-t0,'flux_order':['B','C','D','A'],
       'point_position_order':point_labels,'point_x':np.asarray(x).tolist(),'point_y':np.asarray(y).tolist(),
       'baseline_fluxes':base_flux.tolist(),'substructure_fluxes':full_flux.tolist(),
       'fractional_flux_changes':(full_flux/base_flux-1).tolist(),'S_all_references':score,
       'magnification_baseline':mu_base.tolist(),'magnification_substructure':mu_full.tolist(),
       'point_source_pixel_max_abs_difference':float(np.max(np.abs(point_full-point_base))),
       'subhalos':len(hc.new_subhalos),'fieldhalos':len(hc.new_fieldhalos),
       'lens_components':len(lm['lens_model_list']),
       'adjustment_fraction':float(hc.mass_adjustment_fraction)}
    (OUT/'wgd_fixed_case.json').write_text(json.dumps(r,indent=2))
    print(json.dumps({'wgd_fixed_case_completed':r}),flush=True)


if __name__=='__main__':
    # run.py schedules this after the common-seed maps, alongside the ABC shards.
    from colossus import settings
    settings.BASE_DIR=str(OUT/'diagnostic_cache')
    profile_accuracy()
    wgd_fixed_case()
