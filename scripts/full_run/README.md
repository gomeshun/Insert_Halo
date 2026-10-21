# Thesis-scale full run

These are portable versions of the drivers used for the 2026-10-21 internal
reproduction run. They preserve the scientific calculations in
`Insert_Halo_v1.00` at main commit
`8810e408cb2edfb621ebfbda7aa68d60a4284045`. They do not fix the scientific
issues identified in that audit or produce a validated SIDM constraint.

Only executable drivers, their input settings, and instructions are versioned.
Numerical outputs, figures, fitted systems, logs, and checkpoints stay local.

## Setup and execution

Use Linux x86_64, the Python and uv versions described in the repository README,
Git, and approximately 4 CPUs / 16 GiB RAM. Allow several hours. From the root
of a checkout of this branch:

```bash
uv sync --locked
uv pip check
uv run --locked python scripts/full_run/run.py --plan
uv run --locked python scripts/full_run/run.py --output-dir work/full_run
```

The default output directory is `work/full_run` in this checkout; `work/` is
ignored by Git. Paths are derived from the checkout, not from the original
Cloud machine. `--plan` only prints the schedule and creates no files.
All numerical-library thread limits are set to 1 before numerical imports.
Matplotlib runs headlessly. Colossus uses its supported `settings.BASE_DIR`;
cache directories are confined to the output directory without changing HOME.

The runner first performs a fresh WGD2038 PSO. It then launches the three ABC
shards and one sequence of map/diagnostic jobs, using at most four concurrent
simulation processes. After all ABC shards finish, it runs the additional PSO
controls and the final analysis. Each stage has a separate log. A worker failure
stops the run and terminates the other active workers.

Inspect progress without importing scientific libraries:

```bash
uv run --locked python scripts/full_run/status.py --output-dir work/full_run
```

## Fixed experiment settings

There are no reduced-physics defaults in this runner. A fresh run performs:

| Calculation | Settings |
| --- | --- |
| Initial WGD2038 fit | Original `perform_PSO`: 10 particles, 100 iterations, 10 workers; seed 20261021; fit cosmology H0=67.4, Om0=0.315, Ob0=0.049 |
| ABC | 200 realizations in shards of 67, 67, 66; log10 input halo masses 6–11; log10 host mass 13.5; lens/source redshifts 0.23/0.78 |
| ABC priors | Log-uniform sigma_0 in [0.001, 100] cm²/g and vel_scale in [1, 1000] km/s; other parameters fixed |
| Flux comparison | Data threshold log10(flux)=0.25; connected-component order B,C,D,A; B is the primary reference; S is also recorded for all four references |
| ABC geometry | Cone width and center from the original WGD image configuration; approximately 5.600922395 arcsec wide |
| Common-seed convergence maps | sigma_0=1,10,100 cm²/g, w=100 km/s; seed 20261021 for each case |
| Independent-seed convergence maps | Same three cases, with seeds 20261022, 20261023, 20261024 respectively |
| Map geometry | lens/source redshifts 0.5/2; log10 input halo masses 6–10; log10 host mass 13; width 6 arcsec |
| Map evaluation | All lens components in the original multi-plane `LensModel.kappa`; 300×300 rays; each coordinate is `linspace(-3,3,300)` arcsec; batches of 30 rows |
| Shared halo settings | SIDM; Einasto spatial model; sigmalogc=0.13; fixed input subhalo mass range; original defaults for remaining parameters |
| Shared SASHIMI settings | dz=0.01, N_ma=500, N_herm=20, zmax=10 |
| Gaussian decomposition | 50 Gaussians, P=8; original production bounds GD_min=1e-3, GD_max=1e7 |
| Profile diagnostics | Full halo populations at sigma_0/w=100/25 and 10/24.33; production and calibration GD_max=1e9 errors at R=0.001 kpc; independent two-sided LOS integral |
| Fixed WGD diagnostic | sigma_0=100, w=24.33, seed 20261021; compares full images, point-source-only images, and magnifications |
| Additional macro controls | Three fresh original PSO fits, seeds 20261022–20261024 |

Input halo mass limits refer to the simulation inputs; stripped subhalo masses
can be smaller. The host is excluded from the substructure lens model and from
the convergence maps. The existing ABC macro model is kept fixed across its
200 realizations. Original mass-sheet, cosmology, amplitude, profile, and
Gaussian-filtering behavior is retained.

For ABC shard `s`, the initial NumPy prior RNG seed is `20261021 + 1000*s`.
Realization `i` receives spatial seed `20261021 + 1000000*s + i`.
The NumPy RNG state after each realization is checkpointed and restored on
resume; the interaction of the original simulation with that RNG is retained.
`IH_ABC.run_ABC(N_sim=1)` and its original resets implement the loop. All 200
scores are combined before the original `sort_sum_stats(acc_fraction=0.2)`
retains the best 40. Shards are not independently accepted.

The thesis did not supply its original random seeds, a complete driver, the
original dependency versions, or numerical output arrays. The primary ABC
reference image, cone width, and exact map pixel centers above are explicit
choices made for this reproduction. Fixed seeds document this experiment;
they do not recover the thesis's unknown realization or guarantee identical
floating-point outputs on other hardware.

## Restart and reanalysis

A fresh run refuses to overwrite a nonempty directory. To continue an
interrupted run made by this runner:

```bash
uv run --locked python scripts/full_run/run.py --resume --output-dir work/full_run
```

`run_manifest.json` records the source/dependency fingerprints, driver hashes,
versions, commands, stage status, and the fitted-system hash. Resume rejects
changed drivers, source, environment, or PSO fit. ABC resumes after the last
atomic result/RNG checkpoint. Completed map cases are reused after checking
their configuration and data file. An interrupted map grid is recomputed from
its saved lens model; incomplete diagnostic or control stages restart.
Completed stages must still have their output files.

To rerun just selection, statistics, and figures on complete local outputs:

```bash
uv run --locked python scripts/full_run/run.py --analyze-only --output-dir work/full_run
```

This also accepts the output layout from the original Cloud run, without
requiring a new runner manifest. It needs all three ABC JSON reports, six map
JSON/NPZ pairs, profile and fixed-WGD diagnostic JSON reports, and the three
baseline-control reports in `baseline_controls.json`. The optional original
`paper_image_statistics.json` is consumed if present; digitizing the attached
thesis images is outside these simulation drivers.

Individual workers can also be invoked with uv. Set
`INSERT_HALO_FULL_RUN_DIR` to the same absolute result directory first and run
them in the order enforced by `run.py`. Worker scripts no longer wait
indefinitely for other jobs; the runner manages these dependencies.

## Files and provenance

| Driver | Local output |
| --- | --- |
| `fullrun_pso.py` | `fitted_wgd.pkl`, `original_pso_result.json` |
| `fullrun_abc.py` | Per-shard JSON, atomic result/RNG pickle checkpoint, baseline image, first-halo arrays |
| `fullrun_maps.py`, `fullrun_independent_maps.py` | Halo arrays, cached lens models, full kappa grids, JSON progress and timing |
| `fullrun_diagnostics.py` | Profile error/normalization report, fixed-WGD flux and magnification report, image arrays |
| `fullrun_baseline_controls.py` | Three additional PSO baselines |
| `fullrun_analysis.py` | Global accepted samples, all 200 scores, summary statistics, ABC and convergence figures |

`source_manifest.json` fingerprints the original 17 research Python modules,
the research notebook, and `uv.lock`. It also records SHA-256 hashes of the
seven original Cloud drivers before portability edits. The portable versions
replace machine-specific paths and external provenance-file dependencies,
centralize cache/thread setup, add orchestration and guarded reuse, and make
PSO errors return a failing exit code. Scientific calls and parameters are
preserved. Runtime provenance records the historical main commit separately
from the commit containing these drivers.

The existing CI still performs a fresh PSO and five reduced ABC trials. Fast
CLI/orchestration checks run on push as well. This several-hour full run is
started manually with the command above.
