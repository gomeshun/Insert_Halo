import json
from pathlib import Path
import pickle
import sys
import time
import traceback
from fullrun_common import OUT, configure
configure()
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from WGD2038_setup_class import WGD2038
from colossus import settings as colossus_settings
colossus_settings.BASE_DIR=str(OUT/'cache')

def main():
    out={};t0=time.perf_counter()
    np.random.seed(20261021)
    system=WGD2038()
    try:
        print('Calling original WGD2038.perform_PSO: 10 particles, 100 iterations, threadCount=10',flush=True)
        system.perform_PSO(FlatLambdaCDM(H0=67.4,Om0=.315,Ob0=.049))
        out['status']='completed'
        with (OUT/'fitted_wgd.pkl').open('wb') as f:pickle.dump(system,f)
    except Exception as exc:
        out['status']='error';out['exception']=repr(exc);out['traceback']=traceback.format_exc()
        out['fitted_parameters_present']=system.fitted_model_kwargs is not None
        if system.fitted_model_kwargs is not None:
            with (OUT/'partial_fitted_wgd_kwargs.pkl').open('wb') as f:
                pickle.dump({'kwargs_model':system.kwargs_model,'fitted_model_kwargs':system.fitted_model_kwargs},f)
    out['seconds']=time.perf_counter()-t0
    (OUT/'original_pso_result.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out),flush=True)
    if out['status'] != 'completed':
        raise RuntimeError(out['traceback'])

if __name__=='__main__':
    main()
