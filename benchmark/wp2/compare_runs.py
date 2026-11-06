#!/usr/bin/env python3
"""Bitwise comparison of two run records (P1/P2/P4 JSON): every number and string except timestamps, wall times,
command line and output path. Used for amendment A9 (the glibc malloc setting must not change any result).

Usage:  python benchmark/wp2/compare_runs.py A.json B.json [C.json D.json ...]   (pairs)
"""
from __future__ import annotations

import json
import sys

SKIP = {'start_utc', 'end_utc', 'wall_s', 'argv', 'args', 'environment', 't_equilibrium_s', 't_s'}


def walk(a, b, path=''):
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k not in SKIP:
                out += walk(a.get(k), b.get(k), f'{path}/{k}')
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [(path, f'length {len(a)} vs {len(b)}')]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += walk(x, y, f'{path}[{i}]')
        return out
    return [] if (a == b and type(a) is type(b)) else [(path, a, b)]


def count(a):
    if isinstance(a, dict):
        return sum(count(v) for k, v in a.items() if k not in SKIP)
    if isinstance(a, list):
        return sum(count(v) for v in a)
    return 1


def main():
    args = sys.argv[1:]
    res = []
    for pa, pb in zip(args[0::2], args[1::2]):
        a, b = json.load(open(pa)), json.load(open(pb))
        d = walk(a, b)
        res.append(dict(a=pa, b=pb, values_compared=count(a), differences=len(d), first=[str(x) for x in d[:5]]))
        print(f'{pa}  vs  {pb}: {count(a)} values compared, {len(d)} differ', d[:3] if d else '')
    return 1 if any(r['differences'] for r in res) else 0


if __name__ == '__main__':
    sys.exit(main())
