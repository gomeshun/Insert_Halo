"""Fast checks for the full-run entry point; these never start a simulation."""

import json
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / 'scripts' / 'full_run' / 'run.py'


def invoke(*args):
    return subprocess.run([sys.executable,str(RUNNER),*map(str,args)],
                          capture_output=True,text=True,timeout=15)


def test_plan_preserves_full_scale_and_does_not_create_outputs(tmp_path):
    out = tmp_path/'new-results'
    result = invoke('--plan','--output-dir',out)
    assert result.returncode == 0,result.stderr
    plan = json.loads(result.stdout)
    assert plan['abc_trials'] == 200
    assert plan['shard_totals'] == [67,67,66]
    assert plan['map_shape'] == [300,300]
    assert len(plan['independent_map_cases']) == 3
    assert not out.exists()


def test_fresh_run_refuses_to_overwrite_existing_results(tmp_path):
    marker = tmp_path/'fitted_wgd.pkl'; marker.write_bytes(b'existing-fit')
    result = invoke('--output-dir',tmp_path)
    assert result.returncode != 0
    assert 'not empty' in result.stderr
    assert marker.read_bytes() == b'existing-fit'
    assert sorted(p.name for p in tmp_path.iterdir()) == ['fitted_wgd.pkl']


def test_resume_requires_a_manifest_before_starting_workers(tmp_path):
    result = invoke('--resume','--output-dir',tmp_path)
    assert result.returncode != 0
    assert 'requires run_manifest.json' in result.stderr
    assert list(tmp_path.iterdir()) == []


def load_runner(monkeypatch):
    monkeypatch.syspath_prepend(str(RUNNER.parent))
    spec = importlib.util.spec_from_file_location('full_run_runner',RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    import fullrun_common
    return module,fullrun_common


def test_resume_rejects_a_replaced_pso_fit(tmp_path,monkeypatch):
    module,common = load_runner(monkeypatch)
    runner = module.Runner(tmp_path,False,common)
    fitted = tmp_path/'fitted_wgd.pkl'; fitted.write_bytes(b'original fit')
    runner.record['fit_sha256'] = common.sha256(fitted); runner.save()
    fitted.write_bytes(b'different fit')
    with pytest.raises(RuntimeError,match='PSO fit changed'):
        module.Runner(tmp_path,True,common)


def test_worker_failure_terminates_other_active_workers(tmp_path,monkeypatch):
    module,common = load_runner(monkeypatch)
    fake_plan = module.plan()
    fake_plan['commands'] = {
        name:[sys.executable,'-c','import time; time.sleep(30)']
        for name in module.STAGES
    }
    fake_plan['commands']['pso'] = [sys.executable,'-c',
        'from pathlib import Path; import json; p=Path('+repr(str(tmp_path))+'); '
        '(p/"fitted_wgd.pkl").write_bytes(b"test fit"); '
        '(p/"original_pso_result.json").write_text(json.dumps({"status":"completed"}))']
    fake_plan['commands']['abc-0'] = [sys.executable,'-c','raise SystemExit(7)']
    monkeypatch.setattr(module,'plan',lambda:fake_plan)
    runner = module.Runner(tmp_path,False,common)
    with pytest.raises(RuntimeError,match='abc-0 failed'):
        runner.run()
    assert runner.record['status'] == 'error'
    assert runner.record['stages']['abc-0']['exit_code'] == 7
    assert all(job['proc'].poll() is not None for job in runner.jobs)
