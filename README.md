# Insert_Halo

Insert_Halo combines the publicly available [SASHIMI-SIDM](https://github.com/shinichiroando/sashimi-si) and [pyHalo](https://github.com/dangilman/pyHalo) codes. Dark-matter subhalos generated with SASHIMI are inserted into pyHalo realizations so that they can be used in gravitational-lensing analyses with lenstronomy.

> [!WARNING]
> Insert_Halo is under active development. The code may still contain bugs and interfaces may change.

## Repository layout and versioning

The repository root contains the current version of the code. Historical versions should be recovered from Git history (and, for future releases, Git tags) rather than by keeping a separate source directory for every version.

The current root layout is based on the former `Insert_Halo_v1.00` snapshot.

## Requirements

- Python 3.12
- Runtime dependencies listed and pinned in `pyproject.toml`

The dependency pins in `pyproject.toml` define the reproducible baseline for the v1.00 code snapshot. Dependency updates should be made deliberately and reviewed through Git rather than silently picking up newer releases.

## Installation

Clone the repository and create a Python 3.12 virtual environment:

```bash
git clone https://github.com/gomeshun/Insert_Halo.git
cd Insert_Halo
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

A simple import check is:

```bash
python -c "from IH_main import Insert_Halo; print('Insert_Halo import OK')"
```

## Minimal usage

```python
from IH_main import Insert_Halo

insert_halo = Insert_Halo()
```

The main workflow is exposed through `Insert_Halo.run_insert_halo`. See its docstring and `IH_testing.ipynb` for the code-level interface while the project documentation is being expanded.

## Main external projects

- [SASHIMI-SIDM](https://github.com/shinichiroando/sashimi-si)
- [pyHalo](https://github.com/dangilman/pyHalo)
- [lenstronomy](https://github.com/lenstronomy/lenstronomy)
