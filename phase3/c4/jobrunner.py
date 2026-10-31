"""Checkpointing, bounded job runner for the C-4 diagnostic/TC-1 jobs (replaces the open-ended
`until grep MARKER` chains).

- A job is SKIPPED if its expected output JSON exists and is complete (parsable, has 'end_utc')
  -> verified checkpoint; completed results are never re-run or overwritten.
- Each job runs with a wall-clock limit (subprocess timeout); on timeout the process is killed,
  its log is kept, and the job is recorded as TIMEOUT.
- Every job end (OK / FAILED exit=N / TIMEOUT / SKIPPED) is appended to <out>/RUNS.txt and the
  live state is written to <out>/STATUS.json (pid, job, start, counts).
- A final line 'RUNNER_END ok=.. failed=.. timeout=.. skipped=..' is ALWAYS written to
  <out>/RUNNER_END.txt and stdout (try/finally), so a monitor never waits for a marker that
  cannot appear.
Usage: python phase3/c4/jobrunner.py --jobs JOBS.json --out DIR [--limit-s 2400]
JOBS.json: [{"script": "...", "args": [...], "expect": "<output json path>"}, ...]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds')


def complete(path):
    try:
        return 'end_utc' in json.load(open(path))
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--limit-s', type=float, default=2400.0)
    a = ap.parse_args()
    jobs = json.load(open(a.jobs))
    os.makedirs(os.path.join(a.out, 'logs'), exist_ok=True)
    env = dict(os.environ, DEVSIM_MATH_LIBS='/lib/x86_64-linux-gnu/libblas.so.3:/lib/x86_64-linux-gnu/liblapack.so.3')
    py = os.environ.get('PY', '/home/claude/venv312/bin/python')
    cnt = dict(ok=0, failed=0, timeout=0, skipped=0)
    status_path = os.path.join(a.out, f'STATUS_{os.path.basename(a.jobs)}.json')
    try:
        for i, j in enumerate(jobs):
            tag = '_'.join(x.strip('-') for x in j['args'])
            if complete(j['expect']):
                cnt['skipped'] += 1; res = 'SKIPPED(checkpoint)'; ec = None; st = en = now()
            else:
                st = now()
                log = open(os.path.join(a.out, 'logs', f'{os.path.basename(j["script"])}_{tag}.log'), 'w')
                p = subprocess.Popen([py, j['script'], *j['args'], '--out', a.out], stdout=log, stderr=subprocess.STDOUT, env=env)
                json.dump(dict(runner_pid=os.getpid(), job_index=i, n_jobs=len(jobs), job=j, child_pid=p.pid,
                               start=st, counts=cnt), open(status_path, 'w'), indent=1)
                try:
                    ec = p.wait(timeout=a.limit_s)
                    ok = ec == 0 and complete(j['expect'])
                    res = 'OK' if ok else f'FAILED(exit={ec},output_complete={complete(j["expect"])})'
                    cnt['ok' if ok else 'failed'] += 1
                except subprocess.TimeoutExpired:
                    p.kill(); p.wait(); ec = 'killed'; res = f'TIMEOUT(>{a.limit_s:.0f}s)'; cnt['timeout'] += 1
                log.close(); en = now()
            line = f"job={i} script={j['script']} args='{' '.join(j['args'])}' result={res} start={st} end={en}"
            with open(os.path.join(a.out, 'RUNS.txt'), 'a') as fh:
                fh.write(line + '\n')
            print('JOB', line, flush=True)
    finally:
        end = f"RUNNER_END jobs={os.path.basename(a.jobs)} ok={cnt['ok']} failed={cnt['failed']} timeout={cnt['timeout']} skipped={cnt['skipped']} at={now()}"
        with open(os.path.join(a.out, 'RUNNER_END.txt'), 'a') as fh:
            fh.write(end + '\n')
        json.dump(dict(runner_pid=os.getpid(), finished=True, counts=cnt, end=now()), open(status_path, 'w'), indent=1)
        print(end, flush=True)
    return 0 if cnt['failed'] == cnt['timeout'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
