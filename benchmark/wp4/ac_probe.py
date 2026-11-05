#!/usr/bin/env python3
"""WP4-B2: DEVSIM's built-in small-signal AC analysis on the resolved-well deck R, in double or extended precision.

Same device, mesh and equilibrium as WP2 P1/P2 (benchmark/wp2/rdeck_common.py, unchanged). At equilibrium (V_D = 0) the
drain source V1 carries acreal = 1 V; devsim.solve(type='ac', frequency=F) gives the complex small-signal circuit current.
G_AC = -Re(V1.I) / 1 V is DEVSIM's own linear-response drain conductance; it is compared with the zero-bias network
G_DD on the same mesh (the WP2 reference path). F = 1e-6 Hz, so the capacitive correction to the real part,
O((omega C / G)^2), is below 1e-12 relative for every state.

--extended sets DEVSIM's extended_solver / extended_model / extended_equation before anything is built. In DEVSIM 2.11.2
the AC system is then assembled in float128 but factorized and solved once in double precision (no refinement).

Output: <out>/AC_<id>_<body>_lev<level>_<double|extended>.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'wp2'))
from rdeck_common import ds  # noqa: E402

FLAGS = ('extended_solver', 'extended_model', 'extended_equation')
VERSION = 'wp4-ac-probe-2026-11-04'
FREQ = 1e-6


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--id', required=True)
    for k in ('na', 'tox', 'L', 'vg'):
        ap.add_argument(f'--{k}', type=float, required=True)
    ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--body', choices=['B0', 'B1'], required=True)
    ap.add_argument('--extended', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    prec = 'extended' if a.extended else 'double'
    info = ds.get_parameter(name='info')
    if a.extended:
        if not info.get('extended_precision'):
            raise SystemExit('this DEVSIM build has no extended precision')
        for p in FLAGS:
            ds.set_parameter(name=p, value=True)
    import numpy as np
    import rdeck_common as RC
    from rdeck_common import D, WP
    from rdeck_probe import env
    rec = dict(version=VERSION, start_utc=now(), args=vars(a), argv=sys.argv, environment=env(), precision=prec,
               devsim_info={k: (list(v) if isinstance(v, tuple) else v) for k, v in info.items()}, frequency_Hz=FREQ)
    t0 = time.time()
    _, geo, sol0, b = RC.build(a.L, a.level, a.tox, a.na, a.vg, a.body)
    rec.update(b)
    xs, ys = map(np.array, D.coords(D.SI))
    n0i, n1i, c, asym, psi, n = WP.network(geo)
    M = RC.terminal_masks(geo, xs, ys)
    T = dict(source=np.where(M['source'])[0], drain=np.where(M['drain'])[0])
    if a.body == 'B1':
        T['body'] = np.where(M['body'] & ~M['source'] & ~M['drain'])[0]
    Lap = RC.laplacian(len(xs), n0i, n1i, c)
    rd = RC.terminal_solve(Lap, n0i, n1i, c, T, 'drain')
    G_net = rd['energy']
    t1 = time.time()
    ds.solve(type='ac', frequency=FREQ)
    re = float(ds.get_circuit_node_value(solution='ssac_real', node='V1.I'))
    im = float(ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
    G_ac = -re
    rec.update(network_same_run=dict(G_DD=G_net, residual=rd['residual']),
               ac=dict(I_real=re, I_imag=im, G_AC=G_ac, rel_vs_network=G_ac / G_net - 1.0, t_ac_s=time.time() - t1),
               parameters_after_run={p: (bool(ds.get_parameter(name=p)) if a.extended else False) for p in FLAGS},
               wall_s=time.time() - t0, end_utc=now())
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f'AC_{a.id}_{a.body}_lev{a.level}_{prec}.json')
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(rec, f, indent=1)
    print(f"{a.id} {a.body} lev{a.level} {prec}: G_net={G_net:.6e} G_AC={G_ac:.6e} rel {G_ac / G_net - 1:+.3e} "
          f"(Im {im:+.2e}) wall {rec['wall_s']:.0f}s", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
