#!/usr/bin/env python3
"""End a WP4 queue session at a job boundary so that every attempt is recorded by the queue itself.

1. Wait until the ledger gains its first new line (a running job finished).
2. SIGSTOP the queue process, so it cannot launch further jobs.
3. SIGTERM any child launched after the original ones (its thread records it with exit -15).
4. Wait until every child has exited, then SIGCONT the queue so its threads write their ledger lines.
5. SIGKILL the queue as soon as all expected ledger lines are present; SIGKILL any child it launched in between.
Reason (WP4 report): background tasks are capped at 2 h, shorter than the WP4 queue.

Usage:  python benchmark/wp4/stop_session.py QUEUE_PID LEDGER
"""
from __future__ import annotations

import os
import signal
import subprocess
import sys
import time


def children(pid):
    out = subprocess.run(['ps', '-o', 'pid=,stat=,args=', '--ppid', str(pid)], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        parts = line.split(None, 2)
        if len(parts) >= 2:
            rows.append((int(parts[0]), parts[1], parts[2] if len(parts) > 2 else ''))
    return rows


def nlines(path):
    return sum(1 for _ in open(path)) if os.path.exists(path) else 0


def stamp():
    return time.strftime('%H:%M:%S', time.gmtime())


def main():
    qpid, ledger = int(sys.argv[1]), sys.argv[2]
    n0 = nlines(ledger)
    original = {p for p, _, _ in children(qpid)}
    print(f'{stamp()} queue {qpid}, ledger {n0} lines, running children {sorted(original)}', flush=True)
    while nlines(ledger) == n0:
        time.sleep(0.02)
    os.kill(qpid, signal.SIGSTOP)
    print(f'{stamp()} ledger grew to {nlines(ledger)}; queue stopped (SIGSTOP)', flush=True)
    new = [p for p, st, _ in children(qpid) if p not in original]
    for p in new:
        os.kill(p, signal.SIGTERM)
        print(f'{stamp()} SIGTERM to newly launched child {p}', flush=True)
    while any(st[0] != 'Z' for p, st, _ in children(qpid)):
        time.sleep(2)
    expected = n0 + len(original) + len(new)
    print(f'{stamp()} all children exited; expecting {expected} ledger lines; SIGCONT', flush=True)
    os.kill(qpid, signal.SIGCONT)
    t0 = time.time()
    while nlines(ledger) < expected and time.time() - t0 < 30:
        time.sleep(0.02)
    os.kill(qpid, signal.SIGKILL)
    late = [p for p, st, _ in children(qpid) if st[0] != 'Z']
    for p in late:
        os.kill(p, signal.SIGKILL)
    print(f'{stamp()} queue killed with {nlines(ledger)} ledger lines; children launched after SIGCONT and killed: {late}', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
