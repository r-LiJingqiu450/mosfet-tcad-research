#!/usr/bin/env python3
"""WP8, one state: why the coupled small-signal solve fails (pre-registration benchmark/specs/WP8_*.md).

In one fresh DEVSIM process, for one deck-R state on one mesh level:
  build    WP2 deck and zero-bias equilibrium (benchmark/wp2/rdeck_common.build), and the network reference G_net on
           the same mesh (benchmark/wp2/rdeck_probe.analyse);
  extract  DEVSIM's coupled Jacobian A at equilibrium and the right-hand side of a drain-source step,
           rhs = (b(V1 = 1 mV) - b(V1 = 0)) / 1 mV, exactly as WP2 P3 did (phase3/c4/t1v2.systems);
  M0       G_E (equilibrated SuperLU) and G_M (UMFPACK) recomputed, compared with the WP2 record of the same state;
  M1       the same double-precision system solved by iterative refinement with double-double residuals
           (benchmark/wp8/ddlinalg.py), inner solver the equilibrated SuperLU factorization; M1b the same with plain
           SuperLU as inner solver;
  M2       sensitivity to rounding of the Jacobian entries: three systems with every nonzero of A multiplied by
           (1 + 2^-52 u), u ~ U[-1, 1] (seeds 20261106 + k), each solved as in M1.
G is -x[V1.I] (the circuit branch current per volt), the output DEVSIM's AC analysis reports.

Usage (repository root):
  python benchmark/wp8/ssmech.py --id AR-1 --na 1e18 --tox 2 --L 1 --vg 0.488174 --level 1 --body B1 --out DIR
"""
from __future__ import annotations

import argparse
import gc
import json
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
import scipy.sparse.linalg as spla  # noqa: E402
from phase3.c4.t1v2 import DELTA, EqSolver, systems, umf_solve  # noqa: E402

VERSION = 'wp8-ssmech-2026-11-06'
EPS = 2.0 ** -52
SEED0 = 20261106
N_PERTURB = 3


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def env():
    import devsim
    return dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                devsim=getattr(devsim, '__version__', None), host=platform.node(),
                DEVSIM_MATH_LIBS=os.environ.get('DEVSIM_MATH_LIBS'))


def brief(out, G_net):
    h = out['history']
    return dict(G=out['G'], G_lo=out['G_lo'], rel_vs_G_net=out['G'] / G_net - 1, converged=out['converged'],
                iterations=out['iterations'], omega_dd_final=h[-1]['omega_dd'],
                history=[dict(it=x['it'], G=x['G'], rel_vs_G_net=x['G'] / G_net - 1, omega_dd=x['omega_dd']) for x in h])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--id', required=True)
    ap.add_argument('--na', type=float, required=True)
    ap.add_argument('--tox', type=float, required=True)
    ap.add_argument('--L', type=float, required=True)
    ap.add_argument('--vg', type=float, required=True)
    ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--body', choices=['B0', 'B1'], required=True)
    ap.add_argument('--wp2-record', default=None, help='WP2 P1 record with the small-signal block, for M0')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f'SS_{a.id}_{a.body}_lev{a.level}.json')
    rec = dict(version=VERSION, start_utc=now(), args=vars(a), argv=sys.argv, environment=env())
    t0 = time.time()

    _, geo, sol0, b = RC.build(a.L, a.level, a.tox, a.na, a.vg, a.body)
    rec['onecol_psi_hat'] = b['onecol']['psi_hat']
    rec['equilibrium'] = b['equilibrium']
    net, _ = RP.analyse(geo, sol0, a.body)
    G_net = net['G_DD']
    rec['G_net'] = G_net
    rec['network'] = {k: net[k] for k in ('G_DD', 'G_SD', 'G_DB', 'G_DD_second', 'G_SD_second') if k in net}

    eq = RC.eqnums_capture()
    A, _, b0, bd = systems()
    iV = int(eq[('circuit', 'V1.I')][0])
    rhs = (bd - b0) / DELTA
    nzr = np.nonzero(rhs)[0]
    rec['system'] = dict(n=int(A.shape[0]), nnz=int(A.nnz), index_V1_I=iV, rhs_nonzero_rows=[int(i) for i in nzr[:20]],
                         n_rhs_nonzero=int(len(nzr)), rhs_nonzero_values=[float(rhs[i]) for i in nzr[:20]],
                         b0_inf=float(np.max(np.abs(b0))), A_absmax=float(np.max(np.abs(A.data))),
                         A_absmin=float(np.min(np.abs(A.data[A.data != 0]))))
    rec['t_extract_s'] = time.time() - t0

    # M0: the two double-precision solves of WP2 P3
    t1 = time.time()
    s = EqSolver(A)
    xE = s.solve(rhs)
    G_E = -float(xE[iV])
    xM = umf_solve(A, rhs)
    G_M = -float(xM[iV])
    del xM
    m0 = dict(G_E=G_E, G_M=G_M, rel_G_E=G_E / G_net - 1, rel_G_M=G_M / G_net - 1, eta_E=float(s.eta(A, xE, rhs)))
    if a.wp2_record:
        w = json.load(open(a.wp2_record))
        ss = w['smallsignal']
        m0.update(wp2_G_E=ss['G_E'], wp2_G_M=ss['G_M'], wp2_G_net=w['network']['G_DD'],
                  same_G_E=G_E == ss['G_E'], same_G_M=G_M == ss['G_M'], same_G_net=G_net == w['network']['G_DD'],
                  rel_vs_wp2_G_E=G_E / ss['G_E'] - 1, rel_vs_wp2_G_M=G_M / ss['G_M'] - 1)
    rec['M0'] = m0
    rec['t_M0_s'] = time.time() - t1

    # M1: exact solve of the same double-precision system by refinement with double-double residuals
    t1 = time.time()
    res = DD.DDResidual(A)
    m1 = DD.refine(A, rhs, s.solve, xE, iV, sign=-1.0, residual=res)
    rec['M1'] = brief(m1, G_net)
    G_IR = m1['G']
    rec['t_M1_s'] = time.time() - t1

    # M2: sensitivity to one-ulp perturbations of every Jacobian entry
    t1 = time.time()
    m2 = []
    for k in range(N_PERTURB):
        rng = np.random.default_rng(SEED0 + k)
        Ap = A.copy()
        Ap.data = A.data * (1.0 + EPS * rng.uniform(-1.0, 1.0, A.nnz))
        out = DD.refine(Ap, rhs, s.solve, xE, iV, sign=-1.0)
        m2.append(dict(seed=SEED0 + k, G=out['G'], rel_vs_G_IR=out['G'] / G_IR - 1, rel_vs_G_net=out['G'] / G_net - 1,
                       converged=out['converged'], iterations=out['iterations']))
        del Ap
        gc.collect()
    rec['M2'] = dict(eps=EPS, samples=m2, spread_rel_G_net=max(abs(x['G'] - G_IR) for x in m2) / abs(G_net),
                     all_converged=all(x['converged'] for x in m2))
    rec['t_M2_s'] = time.time() - t1

    # M1b: the same refinement with plain SuperLU (no equilibration) as inner solver
    t1 = time.time()
    del s
    gc.collect()
    try:
        lu = spla.splu(A.tocsc())
        m1b = DD.refine(A, rhs, lu.solve, lu.solve(rhs), iV, sign=-1.0, residual=res)
        rec['M1b'] = brief(m1b, G_net)
        rec['M1b']['rel_vs_G_IR'] = m1b['G'] / G_IR - 1
        del lu
    except Exception as e:  # recorded, not fatal
        rec['M1b'] = dict(error=f'{type(e).__name__}: {e}')
    rec['t_M1b_s'] = time.time() - t1

    rec['wall_s'] = time.time() - t0
    rec['end_utc'] = now()
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(rec, f, indent=1)
    m1r = rec['M1']
    print(f"{a.id} lev{a.level} n={rec['system']['n']} G_net={G_net:.10e} G_E rel {m0['rel_G_E']:+.3e} G_M rel {m0['rel_G_M']:+.3e} | "
          f"M1 rel {m1r['rel_vs_G_net']:+.3e} conv={m1r['converged']} it={m1r['iterations']} | M2 spread {rec['M2']['spread_rel_G_net']:.2e} | "
          f"M1b {rec['M1b'].get('rel_vs_G_net', rec['M1b'].get('error'))} | wall {rec['wall_s']:.0f}s", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
