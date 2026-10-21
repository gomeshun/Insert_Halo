"""Run the fixed thesis-scale experiment; orchestration uses only the stdlib."""

import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT_ENV = "INSERT_HALO_FULL_RUN_DIR"
SHARD_TOTALS = (67, 67, 66)
STAGES = {
    "pso": ("fullrun_pso.py",),
    **{f"abc-{i}": ("fullrun_abc.py", "--shard", str(i), "--total", str(n))
       for i, n in enumerate(SHARD_TOTALS)},
    "maps": ("fullrun_maps.py",),
    "diagnostics": ("fullrun_diagnostics.py",),
    "independent-maps": ("fullrun_independent_maps.py",),
    "baseline-controls": ("fullrun_baseline_controls.py",),
    "analysis": ("fullrun_analysis.py",),
}


def plan():
    return {
        "abc_trials": sum(SHARD_TOTALS), "shard_totals": list(SHARD_TOTALS),
        "acceptance_fraction": 0.2, "map_shape": [300, 300],
        "common_map_seed": 20261021,
        "independent_map_cases": [[1,20261022],[10,20261023],[100,20261024]],
        "schedule": ["pso", "abc-0/abc-1/abc-2 in parallel with maps -> diagnostics -> independent-maps",
                     "wait for all ABC shards", "baseline-controls", "analysis"],
        "commands": {name: [sys.executable, str(HERE / args[0]), *args[1:]]
                     for name, args in STAGES.items()},
    }


def utc():
    return time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())


def check_destination(out, resume, analyze_only):
    if analyze_only:
        if not out.is_dir():
            raise RuntimeError("--analyze-only requires an existing result directory")
    elif resume:
        if not (out / "run_manifest.json").is_file():
            raise RuntimeError("--resume requires run_manifest.json from this runner")
    elif out.exists() and any(out.iterdir()):
        raise RuntimeError("Output directory is not empty; use --resume or choose a new directory")


class Runner:
    def __init__(self, out, resume, common):
        self.out, self.common, self.jobs = out, common, []
        self.path = out / "run_manifest.json"
        source = common.verify_source()
        context = {
            "baseline_commit": source["baseline_commit"],
            "research_files_sha256": source["research_files_sha256"],
            "uv_lock_sha256": source["uv_lock_sha256"],
            "driver_sha256": {p.name:common.sha256(p) for p in sorted(HERE.glob("*.py"))},
            "versions": common.versions(),
            "python": sys.version,
            "plan": {k:v for k,v in plan().items() if k != "commands"},
        }
        if resume:
            self.record = json.loads(self.path.read_text())
            if self.record["context"] != context:
                raise RuntimeError("Source, drivers, environment, or experiment plan changed since this run")
            if "fit_sha256" in self.record:
                if common.sha256(out / "fitted_wgd.pkl") != self.record["fit_sha256"]:
                    raise RuntimeError("The PSO fit changed; it cannot be mixed with existing ABC checkpoints")
        else:
            self.record = {"context":context, "execution_commit":common.execution_commit(),
                           "started_utc":utc(), "stages":{}}
        self.record["status"] = "running"
        self.save()

    def save(self):
        self.common.write_json(self.path, self.record)

    def outputs_complete(self, name):
        def report(relative):
            path = self.out / relative
            return json.loads(path.read_text()) if path.is_file() else {}
        if name == "pso":
            return report("original_pso_result.json").get("status") == "completed" and (self.out/"fitted_wgd.pkl").is_file()
        if name.startswith("abc-"):
            i = int(name[-1]); r = report(f"abc_shard_{i}/abc_full.json")
            return r.get("status") == "completed" and len(r.get("trials",[])) == SHARD_TOTALS[i]
        if name in ("maps", "independent-maps"):
            cases = [(1,20261021),(10,20261021),(100,20261021)] if name == "maps" else [(1,20261022),(10,20261023),(100,20261024)]
            return all(report(f"sigma{s}_w100_seed{seed}.json").get("status") == "completed"
                       and (self.out/f"sigma{s}_w100_seed{seed}_kappa.npz").is_file() for s,seed in cases)
        if name == "diagnostics":
            return all((self.out/f"{n}.json").is_file() for n in ("profile_accuracy_full", "wgd_fixed_case"))
        if name == "baseline-controls":
            r = report("baseline_controls.json")
            return isinstance(r,list) and [v["seed"] for v in r] == [20261022,20261023,20261024]
        r = report("fullrun_summary.json")
        return r.get("abc_completed") == 200 and r.get("accepted") == 40

    def start(self, name):
        prior = self.record["stages"].get(name,{})
        if prior.get("status") == "completed":
            if not self.outputs_complete(name):
                raise RuntimeError(f"Completed stage is missing outputs: {name}")
            print(f"Reusing completed stage: {name}", flush=True)
            return None
        command = plan()["commands"][name]
        log_path = self.out / f"{name}.log"
        handle = log_path.open("a", encoding="utf-8")
        handle.write(f"\n{utc()} Starting {json.dumps(command)}\n"); handle.flush()
        proc = subprocess.Popen(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT,
                                start_new_session=True)
        job = {"name":name,"proc":proc,"handle":handle,"done":False}
        self.jobs.append(job)
        self.record["stages"][name] = {"status":"running", "started_utc":utc(),
            "command":command, "log":log_path.name, "attempt":prior.get("attempt",0)+1}
        self.save(); print(f"Started {name}: {log_path}", flush=True)
        return job

    def poll(self):
        for job in self.jobs:
            if job["done"]:
                continue
            code = job["proc"].poll()
            if code is None:
                continue
            job["handle"].close(); job["done"] = True
            name = job["name"]; ok = code == 0 and self.outputs_complete(name)
            self.record["stages"][name].update(status="completed" if ok else "error",
                                               finished_utc=utc(), exit_code=code)
            self.save()
            if not ok:
                raise RuntimeError(f"Stage {name} failed or has incomplete outputs; see {self.out / (name+'.log')}")
            print(f"Completed {name}", flush=True)

    def wait(self, job):
        if job is None:
            return
        while not job["done"]:
            self.poll()  # Also detect a failure in any concurrent ABC shard.
            if not job["done"]:
                time.sleep(0.5)

    def stop(self):
        for job in self.jobs:
            if job["proc"].poll() is None:
                try:
                    os.killpg(job["proc"].pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    job["proc"].wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(job["proc"].pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    job["proc"].wait()
                self.record["stages"][job["name"]]["status"] = "interrupted"
            job["handle"].close()

    def run(self):
        try:
            self.wait(self.start("pso"))
            self.record["fit_sha256"] = self.common.sha256(self.out/"fitted_wgd.pkl")
            self.save()
            abc = [self.start(f"abc-{i}") for i in range(3)]
            for name in ("maps", "diagnostics", "independent-maps"):
                self.wait(self.start(name))
            for job in abc:
                self.wait(job)
            for name in ("baseline-controls", "analysis"):
                self.wait(self.start(name))
            self.common.verify_source()
            self.record.update(status="completed", finished_utc=utc())
            self.save()
        except BaseException:
            self.stop()
            self.record.update(status="error", finished_utc=utc())
            self.save()
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"work"/"full_run")
    parser.add_argument("--plan", action="store_true", help="Print the fixed experiment plan without running or writing files")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true", help="Resume checkpoints made by this runner")
    mode.add_argument("--analyze-only", action="store_true", help="Reanalyze a complete existing output directory")
    args = parser.parse_args()
    if args.plan:
        print(json.dumps(plan(), indent=2)); return
    out = args.output_dir.expanduser().resolve()
    check_destination(out, args.resume, args.analyze_only)
    os.environ[OUTPUT_ENV] = str(out)
    import fullrun_common as common
    common.verify_source()
    common.configure()
    if args.analyze_only:
        with (out/"analysis.log").open("a", encoding="utf-8") as handle:
            subprocess.run(plan()["commands"]["analysis"], cwd=ROOT, stdout=handle,
                           stderr=subprocess.STDOUT, check=True)
        print(f"Analysis completed: {out/'fullrun_summary.json'}")
    else:
        Runner(out, args.resume, common).run()


if __name__ == "__main__":
    main()
