import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import sys
import time
import traceback

from fullrun_common import ROOT, OUT, SRC, BASELINE_COMMIT, SHARD_TOTALS, configure
configure()


def write_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def main():
    global OUT
    import numpy as np
    parser = argparse.ArgumentParser()
    parser.add_argument('--total', type=int)
    parser.add_argument('--shard', type=int, choices=range(3), default=0)
    args = parser.parse_args()
    if args.total is None:
        args.total = SHARD_TOTALS[args.shard]
    if not 1 <= args.total <= SHARD_TOTALS[args.shard]:
        parser.error('--total must be between 1 and the shard target')
    base_out = OUT
    OUT = base_out / f'abc_shard_{args.shard}'
    OUT.mkdir(parents=True, exist_ok=True)
    from colossus import settings
    settings.BASE_DIR = str(OUT / 'cache')
    from IH_ABC_strong_lensing import IH_ABC, priorc
    outpath = OUT / 'abc_full.json'
    with (base_out / 'fitted_wgd.pkl').open('rb') as f:
        system = pickle.load(f)
    system.setup_data_modelplot()
    obs, labels, _, _ = system.get_flux_image_from_data(flux_threshold_log=.25)
    system.setup_lightmodels(add_lens_light=True, add_shapelets=True)
    system.setup_lensmodel(data_only=True)
    system.setup_image_model()
    kw = system.fitted_model_kwargs
    baseline_image = system.image_model.image(kwargs_lens=kw['kwargs_lens'],
        kwargs_source=kw['kwargs_source'], kwargs_lens_light=kw['kwargs_lens_light'], kwargs_ps=kw['kwargs_ps'])
    baseline_flux = system.get_flux_image_from_sim(baseline_image, labels)
    np.savez_compressed(OUT / 'wgd_baseline.npz', image=baseline_image, labels=labels,
                       observed_image=system.multi_band_list[0][0]['image_data'])
    abc = IH_ABC()
    scores = lambda flux: [float(np.linalg.norm(np.asarray(abc.compute_flux_ratios(obs, i)) -
                                               np.asarray(abc.compute_flux_ratios(flux, i)))) for i in range(4)]
    system.reset_image_model()
    config = dict(zlens=.23, zsource=.78, DM_type='SIDM',
        cone_angle_arcsec=float(system.image_width_arcsec), center_x=float(system.x_center),
        center_y=float(system.y_center), precomputed_host=None, Mhost=13.5,
        msub_min=6, msub_max=11, seed=20261021, spatial='einasto', ngaussian=50,
        sh_dz=.01, sh_N_ma=500, sh_N_herm=20, sh_zmax=10., sigmalogc=.13,
        fix_mass_range_sh=True)
    initial = {'status':'running', 'main_commit':BASELINE_COMMIT,
        'configuration':config, 'intended_trials':200, 'reference_image':'B',
        'flux_order':['B','C','D','A'], 'obs_fluxes':obs,
        'baseline_fluxes':baseline_flux, 'baseline_S_all_references':scores(baseline_flux),
        'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in SRC.glob('*.py')},
        'shard':args.shard, 'prior_rng_seed':20261021 + 1000*args.shard,
        'spatial_seed_policy':'20261021 + 1000000 * shard + trial_index', 'trials':[]}
    checkpoint = OUT / 'checkpoint.pkl'
    saved = None
    if checkpoint.exists():
        with checkpoint.open('rb') as f:
            saved = pickle.load(f)
    result = saved['result'] if saved else (json.loads(outpath.read_text()) if outpath.exists() else initial)
    if len(result['trials']) > args.total:
        raise RuntimeError('Checkpoint contains more trials than requested')
    assert result['configuration'] == config
    assert result['source_sha256'] == initial['source_sha256']
    np.random.seed(20261021 + 1000*args.shard)
    if result['trials']:
        abc.summary_statistics = np.array([r['S_all_references'][0] for r in result['trials']])
        abc.drawn_parameter_values = np.array([r['parameters'] for r in result['trials']], dtype=object)
        abc.LM = [r['lens_components'] for r in result['trials']]
        if saved:
            np.random.set_state(saved['rng'])
        else:
            with (OUT/'abc_rng.pkl').open('rb') as f:
                np.random.set_state(pickle.load(f))
    result['shard'] = args.shard
    result['requested_shard_trials'] = args.total
    result['prior_rng_seed'] = 20261021 + 1000*args.shard
    result['spatial_seed_policy'] = initial['spatial_seed_policy']
    result['status'] = 'running'
    write_json(outpath, result)
    print(json.dumps({'full_configuration':config,'baseline_S_all_references':scores(baseline_flux),
                      'continuing_from':len(result['trials'])}), flush=True)
    priors = {'sigma_0':priorc('sigma_0','loguniform',{'min':1e-3,'max':100.}),
              'vel_scale':priorc('vel_scale','loguniform',{'min':1.,'max':1000.})}
    try:
        for i in range(len(result['trials']), args.total):
            started = time.perf_counter()
            trial_config = dict(config, seed=20261021 + 1000000*args.shard + i)
            abc.run_ABC(prior_dict=priors, tolerance=None, system_class=system,
                obs_fluxes=obs, flux_positions=labels, flux_index_ref=0,
                acceptance_fraction=.2, N_sim=1, **trial_config)
            hc = abc.IH.base_run
            row = {'trial':i, 'spatial_seed':trial_config['seed'], 'parameters':dict(abc.current_draws), 'seconds':time.perf_counter()-started,
                'S_all_references':scores(abc.sim_flux_values), 'fluxes':list(abc.sim_flux_values),
                'subhalos':len(hc.new_subhalos), 'fieldhalos':len(hc.new_fieldhalos),
                'lens_components':abc.LM[-1],
                'truncation_adjustment_fraction':float(hc.truncation_adjustment_fraction),
                'mass_adjustment_fraction':float(hc.mass_adjustment_fraction),
                'max_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
            assert np.isfinite(row['S_all_references']).all()
            result['trials'].append(row)
            result['total_trial_seconds'] = sum(t['seconds'] for t in result['trials'])
            if len(result['trials']) >= 5:
                accepted, values, draws = abc.sort_sum_stats(acc_fraction=.2)
                result['accepted_indices'] = accepted.tolist()
                result['accepted_scores'] = values.tolist()
            with (OUT/'checkpoint.pkl.tmp').open('wb') as f:
                pickle.dump({'result':result,'rng':np.random.get_state()}, f)
            (OUT/'checkpoint.pkl.tmp').replace(checkpoint)
            write_json(outpath, result)
            print(json.dumps({'completed_trial':row}), flush=True)
            if i == 0:
                arrays = {}
                for name, halos in [('sub',hc.new_subhalos),('field',hc.new_fieldhalos)]:
                    for key in ['mass','z','x','y','t_over_tc']:
                        arrays[f'{name}_{key}'] = np.asarray([getattr(h,key,np.nan) for h in halos])
                    for key in ['rs','rc','rt','rhos']:
                        arrays[f'{name}_{key}'] = np.asarray([h.profile_args.get(key,np.nan) for h in halos])
                np.savez_compressed(OUT/'abc_first_halos.npz',**arrays)
            # These are the same resets used between iterations of run_ABC.
            abc.reset_IH()
            abc.reset_sim_fluxes()
            del hc
            gc.collect()
        result['status'] = 'completed'
        write_json(outpath,result)
    except Exception:
        result['status']='error'
        result['traceback']=traceback.format_exc()
        write_json(outpath,result)
        raise


if __name__ == '__main__':
    main()
