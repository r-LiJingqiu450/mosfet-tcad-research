#!/usr/bin/env python3
"""WP2 production run queue (pre-registration section 10).

Builds the pre-registered job list (A-R and held-out states, paths P1/P3, P2, P4), runs at most two jobs at a time
under a memory budget, and appends one ledger line per run (command, start/end, exit code, wall time, peak RSS) to
benchmark/wp2/results/ledger.jsonl. Jobs with a successful ledger entry and an existing output file are skipped, so
the queue can be relaunched until complete. Failed runs stay in the ledger and keep their logs.

Usage (repository root):  python benchmark/wp2/run_queue.py [--budget-min 100] [--only REGEX] [--list]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, 'results', 'production')
LOGS = os.path.join(HERE, 'logs', 'production')
LEDGER = os.path.join(HERE, 'results', 'ledger.jsonl')
PY = sys.executable
MEM_BUDGET_GB = 5.4   # session 2 on: the bash-tool memory cgroup is 5.83 GiB (kernel OOM report, 2026-11-05 01:59 UTC); was 7.0
WORKERS = 2

DEV = {'dev1': dict(na='1e18', tox='2'), 'dev2': dict(na='1e16', tox='20')}
AR = [('AR-1', 'dev1', '1', '0.488174'), ('AR-2', 'dev1', '1', '0.621433'), ('AR-3', 'dev1', '1', '0.7835'),
      ('AR-4', 'dev2', '2', '0.488174'), ('AR-5', 'dev2', '2', '0.621433'), ('AR-6', 'dev2', '2', '0.7835')]
HO = [('HO-1', 'dev1', '1', '0.588534'), ('HO-2', 'dev1', '1', '0.686557'), ('HO-3', 'dev2', '2', '0.588534'),
      ('HO-4', 'dev2', '2', '0.686557'), ('HO-5', 'dev1', '2', '0.488174')]

# peak-memory estimates (GB) from trial runs: dev2 level 2 4.6 GB, dev2 level 1 with P3 3.0 GB
MEM = {('dev2', 0): 1.2, ('dev2', 1): 3.3, ('dev2', 2): 5.2, ('dev1', 0): 0.6, ('dev1', 1): 1.3, ('dev1', 2): 2.2}


def mem_est(path, dev, L, level=None, mesh=None):
    if path == 'P4':
        big = dev == 'dev2' or L == '2'
        return (3.0 if big else 1.8) if mesh == 'refined' else (1.2 if big else 0.8)
    if path == 'P2' and dev == 'dev2' and level == 2:
        return 5.9    # session 1 attempt reached 5.57 GB before the OOM kill: always run alone
    m = MEM[(dev, level)]
    if dev == 'dev1' and L == '2':
        m *= 1.4
    return m


def jobs():
    J = []
    for sid, dev, L, vg in AR:
        base = ['--id', sid, '--na', DEV[dev]['na'], '--tox', DEV[dev]['tox'], '--L', L, '--vg', vg]
        for body in ('B0', 'B1'):
            for lev in (0, 1, 2):
                cmd = [PY, 'benchmark/wp2/rdeck_probe.py', *base, '--level', str(lev), '--body', body, '--out', OUT]
                if lev == 1:
                    cmd.append('--smallsignal')
                J.append(dict(job=f'P1_{sid}_{body}_lev{lev}', cmd=cmd, out=f'P1_{sid}_{body}_lev{lev}.json',
                              mem=mem_est('P1', dev, L, lev)))
        for body, lev in (('B1', 1), ('B1', 2), ('B0', 1)):
            J.append(dict(job=f'P2_{sid}_{body}_lev{lev}', out=f'P2_{sid}_{body}_lev{lev}.json', mem=mem_est('P2', dev, L, lev),
                          cmd=[PY, 'benchmark/wp2/rdeck_dd.py', *base, '--level', str(lev), '--body', body, '--out', OUT]))
        for mesh in ('refined', 'base'):
            J.append(dict(job=f'P4_{sid}_{mesh}', out=f'P4_{sid}_{mesh}.json', mem=mem_est('P4', dev, L, mesh=mesh),
                          cmd=[PY, '-I', 'benchmark/wp2/fv_rdeck.py', *base, '--mesh', mesh, '--out', OUT]))
    for sid, dev, L, vg in HO:
        base = ['--id', sid, '--na', DEV[dev]['na'], '--tox', DEV[dev]['tox'], '--L', L, '--vg', vg]
        for body in ('B0', 'B1'):
            J.append(dict(job=f'P1_{sid}_{body}_lev1', out=f'P1_{sid}_{body}_lev1.json', mem=mem_est('P1', dev, L, 1),
                          cmd=[PY, 'benchmark/wp2/rdeck_probe.py', *base, '--level', '1', '--body', body, '--smallsignal', '--out', OUT]))
        J.append(dict(job=f'P2_{sid}_B1_lev1', out=f'P2_{sid}_B1_lev1.json', mem=mem_est('P2', dev, L, 1),
                      cmd=[PY, 'benchmark/wp2/rdeck_dd.py', *base, '--level', '1', '--body', 'B1', '--out', OUT]))
    return J


def done_jobs():
    ok = set()
    if os.path.exists(LEDGER):
        for line in open(LEDGER):
            r = json.loads(line)
            if r['exit'] == 0 and os.path.exists(os.path.join(OUT, r['out'])):
                ok.add(r['job'])
    return ok


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


LOCK = threading.Lock()


def run_job(j, env):
    log = os.path.join(LOGS, j['job'] + '.log')
    if os.path.exists(log):                        # keep the log of every earlier attempt (failed runs are kept)
        k = 1
        while os.path.exists(os.path.join(LOGS, f"{j['job']}.attempt{k}.log")):
            k += 1
        os.replace(log, os.path.join(LOGS, f"{j['job']}.attempt{k}.log"))
    start, t0 = now(), time.time()
    with open(log, 'w') as f:
        f.write('$ ' + ' '.join(j['cmd']) + '\n'); f.flush()
        p = subprocess.Popen(j['cmd'], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, env=env)
        _, status, ru = os.wait4(p.pid, 0)
        p.returncode = os.waitstatus_to_exitcode(status)
        f.write(f'exit={p.returncode}\n')
    rec = dict(job=j['job'], out=j['out'], cmd=' '.join(c.replace(ROOT + '/', '') for c in j['cmd']), start_utc=start,
               end_utc=now(), exit=p.returncode, wall_s=round(time.time() - t0, 1), maxrss_MB=round(ru.ru_maxrss / 1024),
               mem_estimate_GB=j['mem'])
    with LOCK:
        with open(LEDGER, 'a') as f:
            f.write(json.dumps(rec) + '\n')
    print(f"[{rec['end_utc']}] {j['job']} exit={rec['exit']} wall={rec['wall_s']}s rss={rec['maxrss_MB']}MB", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--budget-min', type=float, default=100.0, help='stop launching new jobs after this many minutes')
    ap.add_argument('--only', default=None)
    ap.add_argument('--list', action='store_true')
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True); os.makedirs(LOGS, exist_ok=True)
    J = jobs()
    if a.only:
        J = [j for j in J if re.search(a.only, j['job'])]
    done = done_jobs()
    todo = sorted([j for j in J if j['job'] not in done], key=lambda j: -j['mem'])
    print(f'{len(J)} jobs, {len(J) - len(todo)} done, {len(todo)} to run', flush=True)
    if a.list:
        for j in todo:
            print(j['job'], j['mem'])
        return 0
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               DEVSIM_MATH_LIBS='/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3')
    running = {}
    t_start = time.time()
    while todo or running:
        for k in [k for k, (th, j) in running.items() if not th.is_alive()]:
            running.pop(k)
        used = sum(j['mem'] for _, j in running.values())
        if todo and len(running) < WORKERS and (time.time() - t_start) / 60 < a.budget_min:
            pick = next((j for j in todo if used + j['mem'] <= MEM_BUDGET_GB), None)
            if pick is None and not running:      # a job above the budget may still run alone
                pick = todo[0]
            if pick is not None:
                todo.remove(pick)
                th = threading.Thread(target=run_job, args=(pick, env), daemon=False)
                th.start(); running[pick['job']] = (th, pick)
                continue
        if not running and (time.time() - t_start) / 60 >= a.budget_min:
            break
        time.sleep(2)
    print(f'queue stop: {len(todo)} jobs left', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
