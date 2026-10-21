"""Paths, cache isolation, and source provenance for the archived full run."""

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "Insert_Halo_v1.00"
OUTPUT_ENV = "INSERT_HALO_FULL_RUN_DIR"
OUT = Path(os.environ.get(OUTPUT_ENV, ROOT / "work" / "full_run")).expanduser().resolve()
SOURCE_MANIFEST = Path(__file__).with_name("source_manifest.json")
BASELINE_COMMIT = "8810e408cb2edfb621ebfbda7aa68d60a4284045"
SHARD_TOTALS = (67, 67, 66)


def configure():
    """Set these before importing numerical libraries, also for PSO children."""
    verify_source()
    cache = OUT / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    for name, value in {
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONPYCACHEPREFIX": str(cache / "pycache"),
        "MPLBACKEND": "Agg",
        "MPLCONFIGDIR": str(cache / "matplotlib"),
        "XDG_CACHE_HOME": str(cache),
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "NUMBA_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    }.items():
        os.environ[name] = value
    sys.dont_write_bytecode = True
    sys.pycache_prefix = str(cache / "pycache")
    sys.path.insert(0, str(SRC))


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_source():
    expected = json.loads(SOURCE_MANIFEST.read_text())
    for relative, digest in expected["research_files_sha256"].items():
        path = ROOT / relative
        if not path.is_file() or sha256(path) != digest:
            raise RuntimeError(f"Research source differs from the archived baseline: {relative}")
    if sha256(ROOT / "uv.lock") != expected["uv_lock_sha256"]:
        raise RuntimeError("uv.lock differs from the archived full-run environment")
    return expected


def execution_commit():
    return subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True
    ).strip()


def versions():
    return {name: importlib.metadata.version(name) for name in (
        "numpy", "scipy", "astropy", "lenstronomy", "pyhalo", "samana",
        "numba", "colossus", "mcfit", "numexpr",
    )}
