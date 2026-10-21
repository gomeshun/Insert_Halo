import argparse
import gc
import json
from pathlib import Path
import pickle
import resource
import sys
import time

from fullrun_common import ROOT, OUT, configure
configure()


def write(path,data):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(path)


def properties(hc,path):
    import numpy as np
    data={}
    for name,halos in [('sub',hc.new_subhalos),('field',hc.new_fieldhalos)]:
        for key in ['mass','z','x','y','t_over_tc','tc','z_acc','r200']:
            data[f'{name}_{key}']=np.asarray([getattr(h,key,np.nan) for h in halos])
        for key in ['rs','rc','rt','rhos']:
            data[f'{name}_{key}']=np.asarray([h.profile_args.get(key,np.nan) for h in halos])
    np.savez_compressed(path,**data)


def run(sigma,w,seed,do_map):
    import numpy as np
    from colossus import settings as colossus_settings
    colossus_settings.BASE_DIR=str(OUT/'map_cache')
    from IH_main import Insert_Halo
    from IH_data_analysis import Insert_Halo_Data_Analysis
    from lenstronomy.LensModel.lens_model import LensModel
    tag=f'sigma{sigma:g}_w{w:g}_seed{seed}'
    jsonpath=OUT/f'{tag}.json'
    lm_path=OUT/f'{tag}_lens.pkl'
    report=json.loads(jsonpath.read_text()) if jsonpath.exists() else {'status':'running'}
    config=dict(zlens=.5,zsource=2.,Mhost=13.,msub_min=6,msub_max=10,cone_angle_arcsec=6.,
        DM_type='SIDM',sigma_0=sigma,vel_scale=w,spatial='einasto',ngaussian=50,seed=seed,
        fix_mass_range_sh=True,sigmalogc=.13,sh_dz=.01,sh_N_ma=500,sh_N_herm=20,sh_zmax=10.)
    if 'configuration' in report and report['configuration'] != config:
        raise RuntimeError(f'Cached map configuration differs: {tag}')
    if report.get('status') == 'completed':
        artifact = OUT / (f'{tag}_kappa.npz' if do_map else f'{tag}_halos.npz')
        if not artifact.is_file():
            raise RuntimeError(f'Completed map is missing its data: {artifact}')
        print(json.dumps({'reusing_completed_case':tag}),flush=True)
        return
    if lm_path.exists() and 'configuration' not in report:
        raise RuntimeError(f'Lens cache has no configuration record: {tag}')
    report['configuration']=config
    if not lm_path.exists():
        np.random.seed(seed);ts=time.perf_counter()
        print(json.dumps({'starting_halo_generation':tag,'configuration':config}),flush=True)
        ih=Insert_Halo();hc=ih.run_insert_halo(**config)
        report['generation_seconds']=time.perf_counter()-ts
        report.update(subhalos=len(hc.new_subhalos),fieldhalos=len(hc.new_fieldhalos),
            truncation_adjustment_fraction=float(hc.truncation_adjustment_fraction),
            mass_adjustment_fraction=float(hc.mass_adjustment_fraction))
        properties(hc,OUT/f'{tag}_halos.npz')
        write(jsonpath,report)
        print(json.dumps({'generation_completed':tag,**report}),flush=True)
        if not do_map:
            report['status']='completed';write(jsonpath,report);return
        ts=time.perf_counter();da=Insert_Halo_Data_Analysis(hc)
        lm=da.prepare_lensmodel(mass_sheet_correction=False,include_subhalos=True,
            include_fieldhalos=True,include_host=False,pre_model_return=True)
        report['lens_conversion_seconds']=time.perf_counter()-ts
        report['lens_components']=len(lm['lens_model_list'])
        with lm_path.open('wb') as f:pickle.dump(lm,f)
        del da,hc,ih
        gc.collect()
        write(jsonpath,report)
        print(json.dumps({'lens_conversion_completed':tag,**report}),flush=True)
    else:
        with lm_path.open('rb') as f:lm=pickle.load(f)
    model=LensModel(lm['lens_model_list'],lens_redshift_list=lm['lens_redshift_array'],
        z_source=2.,multi_plane=True,cosmo=lm['cosmo'])
    grid=np.linspace(-3.,3.,300)
    kappa=np.empty((300,300))
    ts=time.perf_counter()
    # Only the pixel batch size changes: every one of the 300x300 rays uses the
    # original multi-plane LensModel and every generated lens component.
    for lo in range(0,300,30):
        xx,yy=np.meshgrid(grid,grid[lo:lo+30])
        kappa[lo:lo+30]=model.kappa(xx.ravel(),yy.ravel(),lm['kwargs_substructure']).reshape(xx.shape)
        report['rows_done']=min(lo+30,300)
        report['grid_seconds']=time.perf_counter()-ts
        write(jsonpath,report)
        print(json.dumps({'map_rows_completed':tag,'rows':report['rows_done'],
                          'seconds':report['grid_seconds']}),flush=True)
    assert np.isfinite(kappa).all()
    np.savez_compressed(OUT/f'{tag}_kappa.npz',x=grid,kappa=kappa)
    report.update(status='completed',grid_shape=[300,300],kappa_min=float(kappa.min()),
        kappa_max=float(kappa.max()),kappa_mean=float(kappa.mean()),
        kappa_percentiles=np.percentile(kappa,[50,90,95,99,99.9]).tolist(),
        grid_seconds=time.perf_counter()-ts,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    write(jsonpath,report)
    print(json.dumps({'map_completed':tag,**report}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--properties-only',action='store_true')
    args=parser.parse_args()
    if args.properties_only:
        run(100.,25.,20261021,False)
    else:
        for cross_section in [1.,10.,100.]:
            run(cross_section,100.,20261021,True)
