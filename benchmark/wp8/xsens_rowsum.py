#!/usr/bin/env python3
"""WP8 exploratory check X4 (pre-registration section 9, amendment 2): the cancellation in the drain-side n-type rows.

In one fresh DEVSIM process, with A, rhs and the refined x and lambda built exactly as xsens.py, over the
electron-continuity rows of the drain-side n-type region (net doping > 0, x > L/2, contact nodes excluded):
  the relative row sum over the potential columns |sum_{j in psi} A_ij| / sum_{j in psi} |A_ij| (the stored entries summed
  exactly, math.fsum), and over the full row; the range of the refined dpsi_j and dn_j/n_j in that region; for the 10 rows
  with the largest s_i = sum_j |lambda_i A_ij x_j|, the equilibrium n and the largest |A_ij|.

Usage (repository root):
  python benchmark/wp8/xsens_rowsum.py --id AR-1 --na 1e18 --tox 2 --L 1 --vg 0.488174 --level 1 --body B1 --out DIR
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path[:0] = [HERE, os.path.join(ROOT, 'benchmark', 'wp2'), ROOT]

import ddlinalg as DD  # noqa: E402
import rdeck_common as RC  # noqa: E402
import rdeck_probe as RP  # noqa: E402
import xsens as XS  # noqa: E402
from phase3.c4.t1v2 import DELTA, EqSolver, systems  # noqa: E402

VERSION = 'wp8-xsens-rowsum-2026-11-06'


def q(v):
    v = np.asarray(v, dtype=float)
    return dict(min=float(v.min()), median=float(np.median(v)), p99=float(np.percentile(v, 99)), max=float(v.max()), n=int(v.size))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for k in ('--id', '--na', '--tox', '--L', '--vg', '--level', '--body', '--out'):
        ap.add_argument(k, required=True)
    a = ap.parse_args()
    t0 = time.time()
    _, geo, sol0, b = RC.build(float(a.L), int(a.level), float(a.tox), float(a.na), float(a.vg), a.body)
    net, _ = RP.analyse(geo, sol0, a.body)
    eq = RC.eqnums_capture()
    A, _, b0, bd = systems()
    iV = int(eq[('circuit', 'V1.I')][0])
    rhs = (bd - b0) / DELTA
    n = A.shape[0]
    D, ds = RC.D, RC.ds
    ipsi = np.asarray(eq[(D.SI, 'Potential')], dtype=np.int64)
    iele = np.asarray(eq[(D.SI, 'Electrons')], dtype=np.int64)
    iox = np.asarray(eq[(D.OX, 'Potential')], dtype=np.int64)
    is_psi = np.zeros(n, dtype=bool)
    is_psi[ipsi] = True
    is_psi[iox] = True
    xs, ys = map(np.array, D.coords(D.SI))
    netdop = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='NetDop'))
    n_eq = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Electrons'))
    masks = RC.terminal_masks(geo, xs, ys)
    contact = masks['source'] | masks['drain'] | masks['body']
    region = (netdop > 0) & (xs > 0.5 * geo.L_cm) & ~contact          # drain-side n-type, contact nodes excluded

    s = EqSolver(A)
    xE = s.solve(rhs)
    e = np.zeros(n)
    e[iV] = 1.0
    xh, _, hx = XS.refine_to(DD.DDResidual(A), rhs, s.solve, xE)
    lh, _, hl = XS.refine_to(DD.DDResidual(A.transpose().tocsr(copy=True)), e, s.solve_T, s.solve_T(e))

    R = A.tocsr()
    rows = iele[region]
    rel_psi, rel_full, s_row, amax = [], [], [], []
    for i in rows:
        lo, hi = R.indptr[i], R.indptr[i + 1]
        cj, v = R.indices[lo:hi], R.data[lo:hi]
        mp = is_psi[cj]
        sp_abs = math.fsum(np.abs(v[mp]).tolist())
        rel_psi.append(abs(math.fsum(v[mp].tolist())) / sp_abs if sp_abs > 0 else 0.0)
        rel_full.append(abs(math.fsum(v.tolist())) / math.fsum(np.abs(v).tolist()))
        s_row.append(math.fsum(np.abs(lh[i] * v * xh[cj]).tolist()))
        amax.append(float(np.max(np.abs(v))))
    s_row = np.array(s_row)
    top = []
    for k in np.argsort(-s_row)[:10]:
        j = int(np.nonzero(region)[0][k])
        top.append(dict(row=int(rows[k]), x_um=float(xs[j]) * 1e4, y_um=float(ys[j]) * 1e4, n_eq_cm3=float(n_eq[j]),
                        netdop_cm3=float(netdop[j]), max_abs_entry=amax[k], rel_rowsum_psi=rel_psi[k], lambda_i=float(lh[rows[k]]),
                        dpsi_i=float(xh[ipsi[j]])))
    rec = dict(version=VERSION, args=vars(a), n=int(n), region_rows=int(region.sum()), x_omega=hx[-1], lambda_omega=hl[-1],
               rel_rowsum_psi=q(rel_psi), rel_rowsum_full=q(rel_full),
               dpsi_region=q(xh[ipsi[region]]), dn_over_n_region=q(xh[iele[region]] / n_eq[region]),
               lambda_region=q(lh[iele[region]]), n_eq_region=q(n_eq[region]), top_rows=top, G_net=net['G_DD'],
               wall_s=time.time() - t0)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, f'X4_{a.id}_{a.body}_lev{a.level}.json'), 'w') as f:
        json.dump(rec, f, indent=1)
    print(f"{a.id} rows {rec['region_rows']} rel rowsum psi median {rec['rel_rowsum_psi']['median']:.1e} p99 {rec['rel_rowsum_psi']['p99']:.1e} "
          f"max {rec['rel_rowsum_psi']['max']:.1e} | full median {rec['rel_rowsum_full']['median']:.1e} | dpsi [{rec['dpsi_region']['min']:.12f}, "
          f"{rec['dpsi_region']['max']:.12f}] | dn/n max {max(abs(rec['dn_over_n_region']['min']), abs(rec['dn_over_n_region']['max'])):.1e} | "
          f"lambda [{rec['lambda_region']['min']:.6f}, {rec['lambda_region']['max']:.6f}] | wall {rec['wall_s']:.0f}s", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
