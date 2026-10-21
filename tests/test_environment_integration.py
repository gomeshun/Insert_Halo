"""Verify the original WGD2038 PSO and five reduced SIDM ABC trials.

This verifies runtime compatibility, not the scientific validity of the model.
No cached fit, external audit environment, or modifications to the research
modules are used. Run with: uv run --locked pytest -m integration --durations=1
"""

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time
import traceback

import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Insert_Halo_v1.00"
CACHE = ROOT / "work" / "cache"
SEED = 20261021

# Apply before importing numerical libraries, including in PSO worker processes.
# Keep caches writable on headless machines without changing the home directory.
for name, value in {
    "MPLBACKEND": "Agg",
    "MPLCONFIGDIR": str(CACHE / "matplotlib"),
    "XDG_CACHE_HOME": str(CACHE),
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "NUMBA_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "4",
    "PYTHONDONTWRITEBYTECODE": "1",
    "PYTHONPYCACHEPREFIX": str(CACHE / "pycache"),
}.items():
    os.environ[name] = value
sys.dont_write_bytecode = True
sys.pycache_prefix = str(CACHE / "pycache")
# The repository also contains older modules at its root. Always select v1.00.
sys.path.insert(0, str(SOURCE))


def source_hashes():
    return {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(SOURCE.glob("*.py"))
    }


def verify(result):
    import numpy as np
    from astropy.cosmology import FlatLambdaCDM
    from colossus import settings as colossus_settings

    # Colossus documents BASE_DIR as the supported cache-directory setting.
    CACHE.mkdir(parents=True, exist_ok=True)
    colossus_settings.BASE_DIR = str(CACHE)

    from WGD2038_setup_class import WGD2038
    from IH_ABC_strong_lensing import IH_ABC, priorc

    result["versions"] = {
        name: importlib.metadata.version(name)
        for name in (
            "astropy", "colossus", "lenstronomy", "mcfit", "numba",
            "numexpr", "numpy", "pyhalo", "samana", "scipy", "tqdm",
        )
    }
    result["configuration"] = {
        "seed": SEED,
        "pso_particles": 10,
        "pso_iterations": 100,
        "abc_trials": 5,
        "acceptance_fraction": 0.2,
        "z_lens": 0.23,
        "z_source": 0.78,
        "log10_host_mass": 13.5,
        "log10_halo_mass_range": [7, 9],
        "cone_angle_arcsec": 4.0,
        "spatial": "einasto",
        "gaussians": 50,
        "sashimi": {"dz": 0.1, "N_ma": 30, "N_herm": 4, "zmax": 3.0},
        "sigma_logc": 0.13,
        "priors": {
            "sigma_0": {"distribution": "loguniform", "min": 1e-3, "max": 100.0},
            "vel_scale": {"distribution": "loguniform", "min": 1.0, "max": 1000.0},
        },
    }
    before = source_hashes()
    result["source_sha256"] = before
    np.random.seed(SEED)
    system = WGD2038()
    started = time.perf_counter()
    print("Running original WGD2038.perform_PSO (10 particles, 100 iterations)", flush=True)
    system.perform_PSO(FlatLambdaCDM(H0=67.4, Om0=0.315, Ob0=0.049))
    result["pso_seconds"] = time.perf_counter() - started

    # This existing method solves the linear brightness amplitudes omitted by
    # FittingSequence.best_fit(), including the shapelet coefficient array.
    system.setup_data_modelplot()
    observed, labels, _, _ = system.get_flux_image_from_data(flux_threshold_log=0.25)
    system.setup_lightmodels(add_lens_light=True, add_shapelets=True)
    system.setup_lensmodel(data_only=True)
    system.setup_image_model()
    fitted = system.fitted_model_kwargs
    baseline = system.image_model.image(
        kwargs_lens=fitted["kwargs_lens"],
        kwargs_source=fitted["kwargs_source"],
        kwargs_lens_light=fitted["kwargs_lens_light"],
        kwargs_ps=fitted["kwargs_ps"],
    )
    baseline_fluxes = system.get_flux_image_from_sim(baseline, labels)
    if not np.isfinite(baseline).all() or not np.all(np.asarray(observed) > 0):
        raise RuntimeError("Non-finite baseline image or invalid observed fluxes")
    result["observed_fluxes"] = np.asarray(observed).tolist()
    result["baseline_fluxes"] = np.asarray(baseline_fluxes).tolist()
    system.reset_image_model()

    priors = {
        "sigma_0": priorc("sigma_0", "loguniform", {"min": 1e-3, "max": 100.0}),
        "vel_scale": priorc("vel_scale", "loguniform", {"min": 1.0, "max": 1000.0}),
    }
    abc = IH_ABC()
    np.random.seed(SEED)
    started = time.perf_counter()
    print("Running five reduced SIDM ABC trials", flush=True)
    abc.run_ABC(
        prior_dict=priors, tolerance=None, system_class=system,
        zlens=0.23, zsource=0.78, DM_type="SIDM", cone_angle_arcsec=4.0,
        center_x=system.x_center, center_y=system.y_center,
        precomputed_host=None, Mhost=13.5,
        obs_fluxes=observed, flux_positions=labels, flux_index_ref=3,
        acceptance_fraction=0.2, N_sim=5, msub_min=7, msub_max=9,
        seed=SEED, spatial="einasto", ngaussian=50,
        sh_dz=0.1, sh_N_ma=30, sh_N_herm=4, sh_zmax=3.0, sigmalogc=0.13,
    )
    result["abc_seconds"] = time.perf_counter() - started
    if len(abc.summary_statistics) != 5 or not np.isfinite(abc.summary_statistics).all():
        raise RuntimeError("Expected five finite ABC summary statistics")
    accepted, scores, draws = abc.sort_sum_stats(acc_fraction=0.2)
    if len(accepted) != 1 or len(draws) != 1 or not np.isfinite(scores).all():
        raise RuntimeError("ABC selection did not retain one finite-scored sample")
    result.update(
        summary_statistics=abc.summary_statistics.tolist(),
        drawn_parameters=abc.drawn_parameter_values.tolist(),
        accepted_indices=accepted.tolist(),
        accepted_scores=scores.tolist(),
        lens_components_per_trial=abc.LM,
    )
    if source_hashes() != before:
        raise RuntimeError("Research source files changed during verification")


@pytest.mark.integration
def test_original_pso_and_reduced_sidm_abc():
    output_dir = ROOT / "work" / "verification"
    output_dir.mkdir(parents=True, exist_ok=True)
    result = {"status": "running", "python": sys.version, "executable": sys.executable}
    started = time.perf_counter()
    try:
        verify(result)
        result["status"] = "completed"
    except Exception:
        result["status"] = "error"
        result["traceback"] = traceback.format_exc()
        raise
    finally:
        result["total_seconds"] = time.perf_counter() - started
        report = output_dir / "result.json"
        report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(f"Verification {result['status']}: {report}", flush=True)
