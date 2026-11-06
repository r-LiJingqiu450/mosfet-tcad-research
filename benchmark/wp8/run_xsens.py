#!/usr/bin/env python3
"""WP8 exploratory analysis X (pre-registration section 8): run benchmark/wp8/xsens.py for the 11 study states.

Each state runs in a fresh process; a ledger line (kind 'xsens') is appended to benchmark/wp8/results/ledger.jsonl.
Usage (repository root):  python benchmark/wp8/run_xsens.py [--only AR-1]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_wp8 as RW  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only')
    a = ap.parse_args()
    outdir = os.path.join(RW.RES, 'exploratory')
    os.makedirs(outdir, exist_ok=True)
    ledger = os.path.join(RW.RES, 'ledger.jsonl')
    rc = 0
    for s in RW.study_states():
        if a.only and s['id'] != a.only:
            continue
        rec = os.path.join(outdir, f"XS_{s['id']}_B1_lev{s['level']}.json")
        if os.path.exists(rec):
            print(f"{s['id']}: record exists, skipped", flush=True)
            continue
        cmd = [RW.PY, os.path.join(HERE, 'xsens.py'), '--id', s['id'], '--na', repr(s['na']), '--tox', repr(s['tox']),
               '--L', repr(s['L']), '--vg', repr(s['vg']), '--level', str(s['level']), '--body', 'B1', '--out', outdir,
               '--ss-record', os.path.join(RW.RES, 'production', f"SS_{s['id']}_B1_lev{s['level']}.json")]
        t0 = time.time()
        start = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        p = subprocess.run(cmd, cwd=RW.ROOT, env=RW.ENV, capture_output=True, text=True)
        line = dict(kind='xsens', id=s['id'], level=s['level'], start_utc=start, wall_s=round(time.time() - t0, 1),
                    exit=p.returncode, argv=cmd[1:], record=os.path.relpath(rec, RW.ROOT) if os.path.exists(rec) else None,
                    record_sha256=RW.sha(rec) if os.path.exists(rec) else None, stdout_tail=p.stdout[-600:],
                    stderr_tail=p.stderr[-1500:])
        with open(ledger, 'a') as f:
            f.write(json.dumps(line) + '\n')
        print(p.stdout.strip()[-700:] or f"{s['id']}: exit {p.returncode} {p.stderr[-600:]}", flush=True)
        rc |= p.returncode
    return rc


if __name__ == '__main__':
    sys.exit(main())
