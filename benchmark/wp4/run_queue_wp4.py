#!/usr/bin/env python3
"""WP4 run queue (pre-registration benchmark/specs/WP4_Preregistration_External_Paths_2026-11-04.md, section 6).

Jobs: B1 = extended-precision P2 (rdeck_dd_ext.py) on AR-1..AR-6, body B1, level 1 (6 jobs);
      B2 = DEVSIM built-in AC analysis (ac_probe.py) on the same states, double and extended (12 jobs).
Two workers under a 5.4 GB budget, single-threaded BLAS, MALLOC_MMAP_THRESHOLD_=131072 (WP2 amendment A9). One ledger
line per attempt in benchmark/wp4/results/ledger.jsonl; earlier logs are kept as <job>.attemptN.log; completed jobs are
skipped on relaunch; a run still going after MAX_WALL_S is stopped (SIGTERM) and recorded.

Usage (repository root):  python benchmark/wp4/run_queue_wp4.py [--list]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
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
MEM_BUDGET_GB = 5.4
WORKERS = 2
MAX_WALL_S = 4 * 3600
ENV_EXTRA = {'MALLOC_MMAP_THRESHOLD_': '131072'}

AR = [('AR-1', 'dev1', '1', '0.488174'), ('AR-2', 'dev1', '1', '0.621433'), ('AR-3', 'dev1', '1', '0.7835'),
      ('AR-4', 'dev2', '2', '0.488174'), ('AR-5', 'dev2', '2', '0.621433'), ('AR-6', 'dev2', '2', '0.7835')]
DEV = {'dev1': dict(na='1e18', tox='2'), 'dev2': dict(na='1e16', tox='20')}
# (memory GB, expected wall s). Session 1 measured 105 s per extended Newton iteration on the dev2 level-1 mesh (about 58
# iterations per P2 run, 17 per equilibrium); dev1 level 1 has 3.2x fewer nodes. Used only for scheduling.
EST = {('P2X', 'dev1'): (1.2, 2000), ('P2X', 'dev2'): (2.6, 6300), ('ACd', 'dev1'): (0.8, 30), ('ACd', 'dev2'): (2.0, 150),
       ('ACx', 'dev1'): (1.2, 600), ('ACx', 'dev2'): (2.6, 1800)}
FIT = 1.1    # a job is launched only if elapsed + FIT * expected wall fits in the session (background tasks stop at 2 h)


def jobs():
    J = []
    for sid, dev, L, vg in AR:
        base = ['--id', sid, '--na', DEV[dev]['na'], '--tox', DEV[dev]['tox'], '--L', L, '--vg', vg, '--level', '1', '--body', 'B1',
                '--out', OUT]
        m, t = EST[('P2X', dev)]
        J.append(dict(job=f'P2X_{sid}_B1_lev1', out=f'P2X_{sid}_B1_lev1.json', mem=m, t=t,
                      cmd=[PY, 'benchmark/wp4/rdeck_dd_ext.py', *base]))
        for prec, flag in (('double', []), ('extended', ['--extended'])):
            m, t = EST[('ACd' if prec == 'double' else 'ACx', dev)]
            J.append(dict(job=f'AC_{sid}_B1_lev1_{prec}', out=f'AC_{sid}_B1_lev1_{prec}.json', mem=m, t=t,
                          cmd=[PY, 'benchmark/wp4/ac_probe.py', *base, *flag]))
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
    if os.path.exists(log):
        k = 1
        while os.path.exists(os.path.join(LOGS, f"{j['job']}.attempt{k}.log")):
            k += 1
        os.replace(log, os.path.join(LOGS, f"{j['job']}.attempt{k}.log"))
    start, t0 = now(), time.time()
    timed_out = False
    with open(log, 'w') as f:
        f.write('$ ' + ' '.join(j['cmd']) + '\n'); f.flush()
        p = subprocess.Popen(j['cmd'], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, env=env)
        while True:
            pid, status, ru = os.wait4(p.pid, os.WNOHANG)
            if pid:
                break
            if time.time() - t0 > MAX_WALL_S and not timed_out:
                p.terminate(); timed_out = True
            time.sleep(1.0)
        code = os.waitstatus_to_exitcode(status)
        f.write(f'exit={code}\n')
    rec = dict(job=j['job'], out=j['out'], cmd=' '.join(c.replace(ROOT + '/', '') for c in j['cmd']), start_utc=start, end_utc=now(),
               exit=code, wall_s=round(time.time() - t0, 1), maxrss_MB=round(ru.ru_maxrss / 1024), mem_estimate_GB=j['mem'],
               env_extra=ENV_EXTRA, timed_out=timed_out)
    with LOCK:
        with open(LEDGER, 'a') as f:
            f.write(json.dumps(rec) + '\n')
    print(f"[{rec['end_utc']}] {j['job']} exit={code} wall={rec['wall_s']}s rss={rec['maxrss_MB']}MB" + (' TIMEOUT' if timed_out else ''),
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--session-min', type=float, default=1e9, help='launch only jobs expected to finish within this many minutes')
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True); os.makedirs(LOGS, exist_ok=True)
    J = jobs()
    done = done_jobs()
    todo = sorted([j for j in J if j['job'] not in done], key=lambda j: -j['t'])     # longest first
    print(f'{len(J)} jobs, {len(J) - len(todo)} done, {len(todo)} to run', flush=True)
    if a.list:
        for j in todo:
            print(j['job'], j['mem'], j['t'])
        return 0
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               DEVSIM_MATH_LIBS='/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3', **ENV_EXTRA)
    running = {}
    t_start = time.time()
    while todo or running:
        for k in [k for k, (th, j) in running.items() if not th.is_alive()]:
            running.pop(k)
        used = sum(j['mem'] for _, j in running.values())
        left = a.session_min * 60 - (time.time() - t_start)
        if not running and not any(FIT * j['t'] <= left for j in todo):
            break
        if todo and len(running) < WORKERS:
            pick = next((j for j in todo if used + j['mem'] <= MEM_BUDGET_GB and FIT * j['t'] <= left), None)
            if pick is not None:
                todo.remove(pick)
                th = threading.Thread(target=run_job, args=(pick, env), daemon=False)
                th.start(); running[pick['job']] = (th, pick)
                continue
        time.sleep(2)
    print(f'queue stop: {len(todo)} jobs left', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
