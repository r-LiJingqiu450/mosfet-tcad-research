#!/usr/bin/env python3
"""WP8 queue: run benchmark/wp8/ssmech.py for the registered states, one fresh process each, sequentially.

Study: the 11 three-terminal states of tiers A-R and HO on their level-1 meshes, body B1 (the states and mesh of WP2
P-SS and WP4 B2). Smoke (--smoke): one non-study state (dev1, V_G = 0.55 V, level 0) that checks the pipeline; its
outcome is not evaluated. A run whose record exists is skipped. Every run appends a line to the ledger.

Usage (repository root):  python benchmark/wp8/run_wp8.py [--smoke] [--only AR-1]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(HERE, 'results')
PY = sys.executable
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MALLOC_MMAP_THRESHOLD_='131072',
           DEVSIM_MATH_LIBS=os.environ.get('DEVSIM_MATH_LIBS',
                                           '/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3'))
SMOKE = dict(id='SMOKE', na=1e18, tox=2.0, L=1.0, vg=0.55, level=0)


def study_states():
    out = []
    for p in sorted(glob.glob(os.path.join(ROOT, 'benchmark', 'cases', '*.json'))):
        c = json.load(open(p))
        if c['tier'] not in ('A-R', 'HO'):
            continue
        d = c['device']
        out.append(dict(id=d['state_id'], na=d['doping']['N_A_cm3'], tox=d['geometry']['t_ox_nm'], L=d['geometry']['L_um'],
                        vg=d['bias']['V_G_V'], level=1))
    return sorted(out, key=lambda s: s['id'])


def sha(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def run(s, outdir, ledger):
    rec = os.path.join(outdir, f"SS_{s['id']}_B1_lev{s['level']}.json")
    if os.path.exists(rec):
        print(f"{s['id']}: record exists, skipped", flush=True)
        return 0
    cmd = [PY, os.path.join(HERE, 'ssmech.py'), '--id', s['id'], '--na', repr(s['na']), '--tox', repr(s['tox']),
           '--L', repr(s['L']), '--vg', repr(s['vg']), '--level', str(s['level']), '--body', 'B1', '--out', outdir]
    wp2 = os.path.join(ROOT, 'benchmark', 'wp2', 'results', 'production', f"P1_{s['id']}_B1_lev{s['level']}.json")
    if os.path.exists(wp2):
        cmd += ['--wp2-record', wp2]
    t0 = time.time()
    start = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    p = subprocess.run(cmd, cwd=ROOT, env=ENV, capture_output=True, text=True)
    line = dict(id=s['id'], level=s['level'], start_utc=start, wall_s=round(time.time() - t0, 1), exit=p.returncode,
                argv=cmd[1:], record=os.path.relpath(rec, ROOT) if os.path.exists(rec) else None,
                record_sha256=sha(rec) if os.path.exists(rec) else None, stdout_tail=p.stdout[-600:], stderr_tail=p.stderr[-1500:])
    with open(ledger, 'a') as f:
        f.write(json.dumps(line) + '\n')
    print(p.stdout.strip()[-400:] or f"{s['id']}: exit {p.returncode} {p.stderr[-400:]}", flush=True)
    return p.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--only')
    a = ap.parse_args()
    if a.smoke:
        outdir, states = os.path.join(RES, 'smoke'), [SMOKE]
    else:
        outdir, states = os.path.join(RES, 'production'), study_states()
    os.makedirs(outdir, exist_ok=True)
    ledger = os.path.join(RES, 'ledger.jsonl')
    rc = 0
    for s in states:
        if a.only and s['id'] != a.only:
            continue
        rc |= run(s, outdir, ledger)
    return rc


if __name__ == '__main__':
    sys.exit(main())
