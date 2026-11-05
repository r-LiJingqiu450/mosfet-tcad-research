#!/usr/bin/env python3
"""Arm-B verification service for one run.

Listens on a Unix socket in the run's workspace. Each request carries one candidate value of G_DD; the reply is
PASS if the value is positive and within --tol (relative) of the reference, FAIL otherwise. At most --max-calls valid
requests are answered; later ones get "NO CHECKS LEFT". The reference lives only in this process: run_agent.py passes it
on standard input from the private task table, so it does not appear in the process list.
Every request is logged (JSON lines).
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--socket', required=True)
    ap.add_argument('--tol', type=float, default=0.01)
    ap.add_argument('--max-calls', type=int, default=3)
    ap.add_argument('--log', required=True)
    a = ap.parse_args()
    a.reference = float(sys.stdin.readline())
    if os.path.exists(a.socket):
        os.unlink(a.socket)
    srv = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    srv.bind(a.socket)
    srv.listen(4)
    used = 0
    with open(a.log, 'a') as log:
        log.write(json.dumps(dict(t=time.time(), event='start', tol=a.tol, max_calls=a.max_calls)) + '\n'); log.flush()
        while True:
            conn, _ = srv.accept()
            with conn:
                conn.settimeout(10)
                try:
                    raw = conn.recv(256).decode(errors='replace').strip()
                except OSError:
                    continue
                rec = dict(t=time.time(), raw=raw)
                try:
                    v = float(raw)
                    ok_num = v == v and abs(v) != float('inf')
                except ValueError:
                    ok_num = False
                if not ok_num:
                    reply = f'ERROR: not a number: {raw[:60]!r} (no check used)'
                elif used >= a.max_calls:
                    reply = f'NO CHECKS LEFT ({a.max_calls} of {a.max_calls} used)'
                else:
                    used += 1
                    passed = v > 0 and abs(v / a.reference - 1.0) <= a.tol
                    reply = f"{'PASS' if passed else 'FAIL'} (check {used} of {a.max_calls})"
                    rec.update(value=v, verdict='PASS' if passed else 'FAIL', check=used)
                rec['reply'] = reply
                log.write(json.dumps(rec) + '\n'); log.flush()
                try:
                    conn.sendall(reply.encode())
                except OSError:
                    pass


if __name__ == '__main__':
    main()
