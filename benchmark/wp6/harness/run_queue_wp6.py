#!/usr/bin/env python3
"""Execute the WP6 schedule for one stage, two runs at a time, inside one background-task session.

A run is launched only if its worst-case duration (the wall limit plus 120 s) fits in what remains of the session
(--session-min; background tasks stop after 2 hours) and the account's usage, as last reported by Claude Code in a run
transcript (rate_limit_event: 5-hour and 7-day utilization), is below --usage-cap in both windows. Every finished run is
audited (audit.py). Runs already recorded (run.json present) are skipped, so the queue can be relaunched until the stage
is complete.

Usage (repository root):  python benchmark/wp6/harness/run_queue_wp6.py --stage 1 --session-min 118 --usage-cap 0.85
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
WP6 = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(WP6))
RUNS = os.path.join(WP6, 'results', 'runs')
sys.path.insert(0, HERE)
from run_agent import WALL_S  # noqa: E402

PY = sys.executable
WORKERS = 2


def done(run_id):
    return os.path.exists(os.path.join(RUNS, run_id, 'run.json'))


def last_utilization():
    """Latest 5-hour and 7-day utilization reported in any run transcript, as max(5h, 7d) (None if none reported)."""
    best = (0.0, None)
    if not os.path.isdir(RUNS):
        return None
    for rid in os.listdir(RUNS):
        p = os.path.join(RUNS, rid, 'transcript.jsonl.gz')
        if not os.path.exists(p):
            continue
        mt = os.path.getmtime(p)
        with gzip.open(p, 'rt') as f:
            for line in f:
                if '"rate_limit_event"' not in line:
                    continue
                try:
                    info = json.loads(line)['rate_limit_info']
                except (json.JSONDecodeError, KeyError):
                    continue
                w = info.get('unifiedWindows') or {}
                vals = [x.get('utilization') for x in w.values() if isinstance(x, dict)] + [info.get('utilization')]
                vals = [v for v in vals if v is not None]
                if vals and mt >= best[0]:
                    best = (mt, max(vals))
    return best[1]


def run_one(run_id, log):
    t0 = time.time()
    p = subprocess.run([PY, os.path.join(HERE, 'run_agent.py'), '--run-id', run_id], cwd=ROOT, capture_output=True, text=True)
    a = subprocess.run([PY, os.path.join(HERE, 'audit.py'), os.path.join(RUNS, run_id)], cwd=ROOT, capture_output=True, text=True)
    with open(log, 'a') as f:
        f.write(f"[{time.strftime('%H:%M:%S', time.gmtime())}] {run_id} run exit={p.returncode} wall={time.time() - t0:.0f}s "
                f"{p.stdout.strip()[-200:]} {p.stderr.strip()[-300:]}\n  audit exit={a.returncode} {a.stdout.strip()[-400:]} {a.stderr.strip()[-300:]}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', required=True, help="1, 2 or smoke")
    ap.add_argument('--session-min', type=float, default=118)
    ap.add_argument('--usage-cap', type=float, default=0.85)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()
    stage = a.stage if a.stage == 'smoke' else int(a.stage)
    sched = [r for r in json.load(open(os.path.join(HERE, 'schedule.json')))['runs'] if r['stage'] == stage]
    todo = [r['run_id'] for r in sched if not done(r['run_id'])]
    log = a.log or os.path.join(WP6, 'logs', f'queue_stage{a.stage}.log')
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, 'a') as f:
        f.write(f"[{time.strftime('%H:%M:%S', time.gmtime())}] session start: stage {a.stage}, {len(todo)} runs to do, "
                f"usage cap {a.usage_cap}, session {a.session_min} min\n")
    t_start, running, stop_reason = time.time(), {}, None
    while todo or running:
        for k in [k for k, th in running.items() if not th.is_alive()]:
            running.pop(k)
        left = a.session_min * 60 - (time.time() - t_start)
        if todo and len(running) < WORKERS:
            u = last_utilization()
            if u is not None and u >= a.usage_cap:
                stop_reason = f'usage {u:.2f} >= cap {a.usage_cap}'
            elif WALL_S + 120 > left:
                stop_reason = 'session time'
            else:
                rid = todo.pop(0)
                th = threading.Thread(target=run_one, args=(rid, log))
                th.start(); running[rid] = th
                with open(log, 'a') as f:
                    f.write(f"[{time.strftime('%H:%M:%S', time.gmtime())}] launch {rid} (usage {u})\n")
                time.sleep(5)
                continue
            if stop_reason and not running:
                break
            if stop_reason:
                todo = []          # launch nothing more; wait for the running ones
        time.sleep(5)
    with open(log, 'a') as f:
        f.write(f"[{time.strftime('%H:%M:%S', time.gmtime())}] queue stop: {stop_reason or 'stage complete'}; "
                f"{sum(1 for r in sched if not done(r['run_id']))} runs left\n")
    print(f"queue stop: {stop_reason or 'stage complete'}", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
