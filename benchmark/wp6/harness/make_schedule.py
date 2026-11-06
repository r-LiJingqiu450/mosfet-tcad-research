#!/usr/bin/env python3
"""Write the WP6 run schedule (benchmark/wp6/harness/schedule.json), deterministic.

Stage 1 = arm A (no checker), stage 2 = arm B (checker): each 11 tasks x 3 repeats = 33 runs, in a random order drawn
with a fixed seed (numpy default_rng). Run ids W6-001..W6-033 (stage 1) and W6-034..W6-066 (stage 2) carry no task or
arm information; they also name the workspaces. Smoke runs S1 (arm A) and S2 (arm B) use the non-benchmark state SMOKE.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TASKS = [f'T{k:02d}' for k in range(1, 12)]
REPEATS = 3
SEEDS = {1: 20261105, 2: 20261106}


def main():
    runs = [dict(run_id='S1', task='SMOKE', arm='A', repeat=0, stage='smoke'),
            dict(run_id='S2', task='SMOKE', arm='B', repeat=0, stage='smoke')]
    k = 0
    for stage, arm in ((1, 'A'), (2, 'B')):
        pairs = [(t, r) for t in TASKS for r in range(1, REPEATS + 1)]
        order = np.random.default_rng(SEEDS[stage]).permutation(len(pairs))
        for i in order:
            k += 1
            t, rep = pairs[int(i)]
            runs.append(dict(run_id=f'W6-{k:03d}', task=t, arm=arm, repeat=rep, stage=stage))
    out = dict(note='WP6 run schedule; runs are executed in this order, two at a time', seeds=SEEDS, runs=runs)
    with open(os.path.join(HERE, 'schedule.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(f'{len(runs)} runs written')
    return 0


if __name__ == '__main__':
    sys.exit(main())
