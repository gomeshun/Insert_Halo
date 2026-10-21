import json
from pathlib import Path
import pickle
import sys
import time

from fullrun_common import OUT, configure
configure()


def main():
    # run.py waits for all ABC shards before starting these PSO fits.
    import numpy as np
    from astropy.cosmology import FlatLambdaCDM
    from colossus import settings
    settings.BASE_DIR=str(OUT/'baseline_controls_cache')
    from WGD2038_setup_class import WGD2038
    from IH_ABC_strong_lensing import IH_ABC
    rows=[]
    for seed in [20261022,20261023,20261024]:
        ts=time.perf_counter();np.random.seed(seed)
        s=WGD2038();s.perform_PSO(FlatLambdaCDM(H0=67.4,Om0=.315,Ob0=.049))
        s.setup_data_modelplot()
        obs,labels,_,_=s.get_flux_image_from_data(flux_threshold_log=.25)
        s.setup_lightmodels(add_lens_light=True,add_shapelets=True)
        s.setup_lensmodel(data_only=True);s.setup_image_model()
        kw=s.fitted_model_kwargs
        image=s.image_model.image(kwargs_lens=kw['kwargs_lens'],kwargs_source=kw['kwargs_source'],
            kwargs_lens_light=kw['kwargs_lens_light'],kwargs_ps=kw['kwargs_ps'])
        flux=s.get_flux_image_from_sim(image,labels);abc=IH_ABC()
        scores=[float(np.linalg.norm(np.asarray(abc.compute_flux_ratios(obs,i))-
                                     np.asarray(abc.compute_flux_ratios(flux,i)))) for i in range(4)]
        row={'seed':seed,'seconds':time.perf_counter()-ts,'baseline_fluxes':flux,
             'baseline_S_all_references':scores}
        rows.append(row)
        (OUT/'baseline_controls.json').write_text(json.dumps(rows,indent=2))
        print(json.dumps({'baseline_control':row}),flush=True)


if __name__=='__main__':main()
