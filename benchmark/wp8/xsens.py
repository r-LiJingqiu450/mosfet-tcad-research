#!/usr/bin/env python3
"""WP8 exploratory analysis X (pre-registration benchmark/specs/WP8_*.md, section 8, amendment 1). Descriptive only.

In one fresh DEVSIM process, for one deck-R state on one mesh level, with A, rhs, G_net built exactly as ssmech.py:
  x, lambda  the solution x = A^-1 rhs and the adjoint lambda = A^-T e_V1.I, both refined with double-double residuals
             (benchmark/wp8/ddlinalg.py) until the componentwise backward error is <= 1e-24 (at most 12 iterations);
             inner solvers: the equilibrated SuperLU factorization and its transpose
  X2         first-order coefficients c_ij = lambda_i A_ij x_j over the nonzeros: a perturbation dA changes G by
             dG1 = sum_ij c_ij dA_ij/A_ij. kappa_A = sum|c|/|G_net|; absolute scales u sum|c| and u (sum c^2)^(1/2)
             (u = 2^-53); the same quantities from the unrefined double-precision x and lambda, and WP2's kappa_G
  X1         for each M2 seed: the perturbed matrix rebuilt exactly as ssmech.py, dA = A_k - A (exact), the first-order
             term dG1 = lambda^T dA x, the second-order term dG2 = -lambda^T dA z with z = A^-1 (dA x) (refined), and the
             measured G_k - G_IR (re-solved exactly as M2, compared bit for bit with the M2 record when given)
  X3         s_i = sum_j |c_ij| per row, aggregated by equation block and by the class of the row's node (contact
             node, n-type drain side / source side, p-type, oxide, circuit); the 10 largest rows; a 64 x 32 map of the
             silicon-node share
G is -x[V1.I], as in ssmech.py.

Usage (repository root):
  python benchmark/wp8/xsens.py --id AR-1 --na 1e18 --tox 2 --L 1 --vg 0.488174 --level 1 --body B1 \
      --ss-record benchmark/wp8/results/production/SS_AR-1_B1_lev1.json --out DIR
"""
from __future__ import annotations

import argparse
import gc
import json
import math
import os
import platform
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path[:0] = [HERE, os.path.join(ROOT, 'benchmark', 'wp2'), ROOT]

import ddlinalg as DD  # noqa: E402
import rdeck_common as RC  # noqa: E402
import rdeck_probe as RP  # noqa: E402
import scipy  # noqa: E402
from phase3.c4.t1v2 import DELTA, EqSolver, systems  # noqa: E402

VERSION = 'wp8-xsens-2026-11-06'
EPS = 2.0 ** -52          # the M2 perturbation scale (ssmech.py)
U = 2.0 ** -53            # unit roundoff
SEED0 = 20261106
N_PERTURB = 3
OMEGA_TOL = 1e-24
MAX_IT = 12
BLOCKS = ('si_potential', 'si_electrons', 'ox_potential', 'circuit')


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def refine_to(res, b, solve, x0, omega_tol=OMEGA_TOL, max_it=MAX_IT):
    """Iterative refinement with double-double residuals until the componentwise backward error <= omega_tol."""
    b = np.asarray(b, dtype=np.float64)
    xh = np.array(x0, dtype=np.float64)
    xl = np.zeros_like(xh)
    hist = []
    for k in range(max_it + 1):
        rh, rl = res(b, xh, xl)
        r = rh + rl
        om = res.componentwise_backward_error(r, b, xh)
        hist.append(om)
        if om <= omega_tol or k == max_it:
            break
        d = solve(r)
        s, t = DD.two_sum(xh, d)
        xh, xl = DD.fast_two_sum(s, t + xl)
    return xh, xl, hist


def fsum(v):
    return math.fsum(np.asarray(v, dtype=np.float64).tolist())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--id', required=True)
    ap.add_argument('--na', type=float, required=True)
    ap.add_argument('--tox', type=float, required=True)
    ap.add_argument('--L', type=float, required=True)
    ap.add_argument('--vg', type=float, required=True)
    ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--body', choices=['B0', 'B1'], required=True)
    ap.add_argument('--ss-record', default=None, help='the ssmech.py record of the same state (M1, M2 values)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f'XS_{a.id}_{a.body}_lev{a.level}.json')
    import devsim
    rec = dict(version=VERSION, start_utc=now(), args=vars(a), argv=sys.argv,
               environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                                devsim=getattr(devsim, '__version__', None), host=platform.node(),
                                DEVSIM_MATH_LIBS=os.environ.get('DEVSIM_MATH_LIBS')))
    t0 = time.time()

    # Build exactly as ssmech.py
    _, geo, sol0, b = RC.build(a.L, a.level, a.tox, a.na, a.vg, a.body)
    net, _ = RP.analyse(geo, sol0, a.body)
    G_net = net['G_DD']
    eq = RC.eqnums_capture()
    A, _, b0, bd = systems()
    iV = int(eq[('circuit', 'V1.I')][0])
    rhs = (bd - b0) / DELTA
    n = A.shape[0]
    rec['G_net'] = G_net
    rec['system'] = dict(n=int(n), nnz=int(A.nnz), index_V1_I=iV, format=A.format)

    # Row map: block and node of every equation number (equation and unknown numbering coincide in DEVSIM)
    D, ds = RC.D, RC.ds
    block = np.full(n, -1, dtype=np.int64)
    node = np.full(n, -1, dtype=np.int64)
    for code, key in enumerate(((D.SI, 'Potential'), (D.SI, 'Electrons'), (D.OX, 'Potential'))):
        idx = np.asarray(eq[key], dtype=np.int64)
        block[idx] = code
        node[idx] = np.arange(len(idx))
    circ = {k[1]: int(v[0]) for k, v in eq.items() if k[0] == 'circuit'}
    for v in circ.values():
        block[v] = 3
    rec['system']['circuit_rows'] = circ
    rec['system']['unassigned_rows'] = int(np.sum(block < 0))
    xs, ys = map(np.array, D.coords(D.SI))
    netdop = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='NetDop'))
    masks = RC.terminal_masks(geo, xs, ys)
    L = geo.L_cm
    si_class = np.where(netdop > 0, np.where(xs > 0.5 * L, 'n_drain_side', 'n_source_side'), 'p_type').astype(object)
    for t in ('source', 'drain', 'body'):
        si_class[masks[t]] = f'contact_{t}'

    def row_class(i):
        bl = block[i]
        if bl in (0, 1):
            return si_class[node[i]]
        return 'oxide' if bl == 2 else 'circuit'

    rec['t_build_s'] = time.time() - t0

    # Solutions: unrefined (as WP2) and refined
    t1 = time.time()
    s = EqSolver(A)
    xE = s.solve(rhs)
    e = np.zeros(n)
    e[iV] = 1.0
    lE = s.solve_T(e)
    res = DD.DDResidual(A)
    AT = A.transpose().tocsr(copy=True)
    resT = DD.DDResidual(AT)
    xh, xl, hx = refine_to(res, rhs, s.solve, xE)
    lh, ll, hl = refine_to(resT, e, s.solve_T, lE)
    G_x = -float(xh[iV] + xl[iV])
    G_l = -float(lh[iV] + ll[iV])
    # M1 exactly as ssmech.py (its G_IR is the reference point of the M2 shifts)
    m1 = DD.refine(A, rhs, s.solve, xE, iV, sign=-1.0, residual=res)
    G_IR = m1['G']
    rec['solutions'] = dict(x_omega_history=hx, lambda_omega_history=hl, G_from_x=G_x, G_from_lambda=G_l,
                            rel_x_vs_lambda=G_x / G_l - 1, G_IR_M1=G_IR, rel_x_vs_M1=G_x / G_IR - 1,
                            rel_G_IR_vs_G_net=G_IR / G_net - 1, G_E=-float(xE[iV]))
    rec['t_solve_s'] = time.time() - t1

    # X2: first-order coefficients
    t1 = time.time()
    rows = A.indices.astype(np.int64)
    cols = np.repeat(np.arange(n, dtype=np.int64), np.diff(A.indptr))
    c = lh[rows] * A.data * xh[cols]
    ac = np.abs(c)
    sum_abs, sum_sq = fsum(ac), fsum(c * c)
    cE = np.abs(lE[rows] * A.data * xE[cols])
    G_E = -float(xE[iV])
    sc = abs(A) @ np.abs(xE) + np.abs(rhs)
    rec['X2'] = dict(kappa_A=sum_abs / abs(G_net), abs_sum_u=U * sum_abs, abs_rms_u=U * math.sqrt(sum_sq),
                     rel_rms_u=U * math.sqrt(sum_sq) / abs(G_net), signed_sum_rel=fsum(c) / G_net,
                     unrefined=dict(kappa_A_vs_G_net=fsum(cE) / abs(G_net), abs_sum_u=U * fsum(cE),
                                    abs_rms_u=U * math.sqrt(fsum(cE * cE)),
                                    kappa_G_wp2=float(np.abs(lE) @ sc) / abs(G_E)))
    rec['t_X2_s'] = time.time() - t1

    # X1: linearity of the response to the M2 perturbations
    t1 = time.time()
    ssr = json.load(open(a.ss_record)) if a.ss_record else None
    x1 = []
    for k in range(N_PERTURB):
        rng = np.random.default_rng(SEED0 + k)
        Ap = A.copy()
        Ap.data = A.data * (1.0 + EPS * rng.uniform(-1.0, 1.0, A.nnz))
        dA = Ap.data - A.data                                   # exact (Sterbenz)
        nzm = A.data != 0
        rho = dA[nzm] / A.data[nzm]
        ulps = dA[nzm] / np.spacing(np.abs(A.data[nzm]))
        dG1 = fsum(lh[rows] * dA * xh[cols])
        w = np.bincount(rows, weights=dA * xh[cols], minlength=n)
        zh, zl, hz = refine_to(res, w, s.solve, s.solve(w))
        dG2 = -fsum(lh[rows] * dA * zh[cols])
        out = DD.refine(Ap, rhs, s.solve, xE, iV, sign=-1.0)
        meas = out['G'] - G_IR
        item = dict(seed=SEED0 + k, rho_mean=float(rho.mean()), rho_rms=float(np.sqrt(np.mean(rho * rho))),
                    ulp_mean=float(ulps.mean()), ulp_rms=float(np.sqrt(np.mean(ulps * ulps))),
                    G_k=out['G'], converged=out['converged'], omega_dd=out['history'][-1]['omega_dd'],
                    measured_shift=meas, dG1=dG1, dG2=dG2, z_omega=hz[-1],
                    rel_dG1_vs_measured=dG1 / meas - 1 if meas else None,
                    rel_dG12_vs_measured=(dG1 + dG2) / meas - 1 if meas else None,
                    shift_rel_G_net=meas / abs(G_net))
        if ssr:
            mk = ssr['M2']['samples'][k]
            item.update(same_as_M2_record=out['G'] == mk['G'] and G_IR == ssr['M1']['G'])
        x1.append(item)
        del Ap, dA, w, zh, zl
        gc.collect()
    rec['X1'] = x1
    rec['t_X1_s'] = time.time() - t1

    # X3: location
    t1 = time.time()
    s_row = np.bincount(rows, weights=ac, minlength=n)
    row_nnz = np.bincount(rows, minlength=n)
    s_col = np.bincount(cols, weights=ac, minlength=n)
    tot = float(s_row.sum())
    classes = np.array([row_class(i) for i in range(n)], dtype=object)
    by = {}
    for bl in range(4):
        for cl in sorted(set(classes[block == bl])):
            m = (block == bl) & (classes == cl)
            by[f'{BLOCKS[bl]}|{cl}'] = dict(rows=int(m.sum()), share=float(s_row[m].sum()) / tot)
    by_class, by_block, by_col_block = {}, {}, {}
    for cl in sorted(set(classes)):
        by_class[cl] = float(s_row[classes == cl].sum()) / tot
    for bl in range(4):
        by_block[BLOCKS[bl]] = float(s_row[block == bl].sum()) / tot
        by_col_block[BLOCKS[bl]] = float(s_col[block == bl].sum()) / tot
    top = []
    for i in np.argsort(-s_row)[:10]:
        i = int(i)
        t = dict(row=i, block=BLOCKS[block[i]], cls=classes[i], share=float(s_row[i]) / tot, lambda_i=float(lh[i]),
                 x_i=float(xh[i]), n_entries=int(row_nnz[i]), max_abs_entry=float(np.max(np.abs(A.data[rows == i]))))
        if block[i] in (0, 1):
            j = node[i]
            t.update(x_um=float(xs[j]) * 1e4, y_um=float(ys[j]) * 1e4, netdop_cm3=float(netdop[j]))
        top.append(t)
    # map of the silicon-node share, 64 x 32 cells
    s_node = np.zeros(len(xs))
    for bl in (0, 1):
        idx = np.asarray(eq[(D.SI, ('Potential', 'Electrons')[bl])], dtype=np.int64)
        s_node += s_row[idx]
    ex = np.linspace(xs.min(), xs.max(), 65)
    ey = np.linspace(ys.min(), ys.max(), 33)
    H, _, _ = np.histogram2d(xs, ys, bins=[ex, ey], weights=s_node / tot)
    rec['X3'] = dict(total_abs=tot, by_block_and_class=by, by_class=by_class, by_row_block=by_block,
                     by_column_block=by_col_block, top_rows=top, silicon_share=float(s_node.sum()) / tot,
                     map=dict(x_edges_um=[float(v) * 1e4 for v in ex], y_edges_um=[float(v) * 1e4 for v in ey],
                              share=[[float(f'{v:.3e}') for v in row] for row in H]),
                     geometry=dict(L_um=L * 1e4))
    rec['t_X3_s'] = time.time() - t1

    rec['wall_s'] = time.time() - t0
    rec['end_utc'] = now()
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(rec, f, indent=1)
    X2 = rec['X2']
    print(f"{a.id} n={n} G_net={G_net:.4e} x/lambda {rec['solutions']['rel_x_vs_lambda']:+.1e} | kappa_A {X2['kappa_A']:.2e} "
          f"abs_rms {X2['abs_rms_u']:.2e} (unrefined {X2['unrefined']['abs_rms_u']:.2e}, WP2 kappa {X2['unrefined']['kappa_G_wp2']:.1e}) | "
          + ' '.join(f"[{x['measured_shift']:+.3e} dG1 {x['dG1']:+.3e} dG2 {x['dG2']:+.1e}]" for x in x1)
          + f" | drain-side n {by_class.get('n_drain_side', 0):.2f} wall {rec['wall_s']:.0f}s", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
