"""Print current progress without importing any scientific modules."""

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'work'/'full_run')
    out = parser.parse_args().output_dir.expanduser().resolve()
    for i,total in enumerate((67,67,66)):
        p = out/f'abc_shard_{i}'/'abc_full.json'
        r = json.loads(p.read_text()) if p.exists() else {}
        print(f"ABC {i}: {len(r.get('trials',[]))}/{total} {r.get('status','not started')}")
    for p in sorted(out.glob('sigma*_w100_seed*.json')):
        r = json.loads(p.read_text())
        print(f"{p.stem}: {r.get('rows_done',0)}/300 rows {r['status']}")
    p = out/'run_manifest.json'
    if p.exists():
        r = json.loads(p.read_text()); print('Overall:',r['status'])
        for name,record in r['stages'].items():
            print(f"  {name}: {record['status']}")


if __name__ == '__main__':
    main()
