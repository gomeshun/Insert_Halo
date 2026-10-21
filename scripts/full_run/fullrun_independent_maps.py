import json
from pathlib import Path
import time

from fullrun_common import OUT, configure
configure()
if __name__=='__main__':
    from fullrun_maps import run
    for sigma,seed in [(1.,20261022),(10.,20261023),(100.,20261024)]:
        run(sigma,100.,seed,True)
