"""POST-HOC DIAGNOSTIC (does not change any T1 verdict; Research_08 section 3.5).

Purpose: explain the single NON-RESOLVING (class I) T1 point (L = 10 um, level 3, R-I(a),
V_G = 2.0 V), where the L1 path (SciPy SuperLU on the exported Newton matrix) differs from
DEVSIM's own UMFPACK Newton step by 2.7e-6 relative. Questions:
  (1) what is the normwise backward error of the SuperLU solution used by L1?
  (2) does iterative refinement (same factorization, residual in double precision) move
      G_L1 to G_NS?
  (3) does row/column equilibration before factorization do the same?
Nothing here is used to re-evaluate T1. Output: JSON with all numbers.
Usage: python phase3/c3/l1_refine_diag.py --L 10 --level 3 --mode a --vg 2.0 --out DIR
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import devsim as ds

from phase3 import runlog
from phase3.c1 import deck as D
from phase3.c1.c1b_jacobian import export
from phase3.c1.corner_study import eqnums
from phase3.c3.t1 import DELTA_L1, Geo, equilibrate, g_gca, newton_step_G


def backward_error(A, x, b):
    """Returns (normwise inf-norm backward error, componentwise Oettli-Prager backward error,
    residual). The normwise value is dominated by the largest rows of this badly scaled
    system (run1 showed ~1e-37) and is kept only for completeness; the componentwise value
    max_i |r_i| / (|A||x| + |b|)_i is invariant under row scaling and is the one interpreted."""
    r = b - A @ x
    nA = abs(A).sum(axis=1).max()
    nw = float(np.abs(r).max() / (nA * np.abs(x).max() + np.abs(b).max()))
    den = abs(A) @ np.abs(x) + np.abs(b)
    m = den > 0
    cw = float(np.max(np.abs(r[m]) / den[m])) if m.any() else 0.0
    return nw, cw, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--mode', default='a'); ap.add_argument('--vg', type=float, default=2.0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    rec = runlog.header(dict(vars(a), note='post-hoc diagnostic; not used for any verdict'))
    geo = Geo(a.L, a.level)
    D.build_mesh_explicit(geo.xs, geo.y_ox, geo.y_si, corner='tie')
    D.build_physics()
    ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
    sol0, rec['equilibrium'] = equilibrate(geo, a.vg, a.mode)
    G, _ = g_gca(geo, sol0); rec['G_GCA_h'] = G
    eq = eqnums(); iV = int(eq[('circuit', 'V1.I')][0])
    A, b0, *_ = export()
    ds.circuit_alter(name='V1', value=DELTA_L1)
    _, bd, *_ = export()
    ds.circuit_alter(name='V1', value=0.)
    rhs = bd - b0
    A = A.tocsc()
    out = dict(n=int(A.shape[0]), nnz=int(A.nnz))
    # B3-consistent metric (phase3/c1/c1b_jacobian.py): normwise backward error of the
    # row+column max-equilibrated system. (Componentwise omega is degenerate = 1 on decoupled
    # contact rows with b = 0; the raw normwise value is dominated by the largest rows.)
    absA = abs(A)
    Dr_ = 1.0 / np.asarray(absA.max(axis=1).todense()).ravel()
    Ar_ = sp.diags(Dr_) @ A
    Dc_ = 1.0 / np.asarray(abs(Ar_).max(axis=0).todense()).ravel()
    Ae_ = (Ar_ @ sp.diags(Dc_)).tocsr(); be_ = Dr_ * rhs; nAe_ = float(abs(Ae_).sum(axis=1).max())

    def eta(x):
        xe = x / Dc_
        return float(np.max(np.abs(Ae_ @ xe - be_)) / (nAe_ * np.max(np.abs(xe)) + np.max(np.abs(be_))))
    # (1)+(2) plain SuperLU as in t1.l1, then iterative refinement
    lu = spla.splu(A)
    x = lu.solve(rhs); steps = []
    for k in range(4):
        nw, cw, r = backward_error(A, x, rhs)
        steps.append(dict(step=k, backward_error_normwise=nw, backward_error_componentwise=cw,
                          backward_error_equilibrated=eta(x),
                          G=-float(x[iV]) / DELTA_L1, r0=-float(x[iV]) / DELTA_L1 / G - 1))
        x = x + lu.solve(r)
    out['superlu_refinement'] = steps
    # (3) equilibrated factorization: Dr A Dc, row then column max-scaling
    dr = 1.0 / np.asarray(abs(A).max(axis=1).todense()).ravel()
    As = sp.diags(dr) @ A
    dc = 1.0 / np.asarray(abs(As).max(axis=0).todense()).ravel()
    As = (As @ sp.diags(dc)).tocsc()
    lus = spla.splu(As)
    y = lus.solve(dr * rhs); xs = dc * y
    nw, cw, _ = backward_error(A, xs, rhs)
    out['superlu_equilibrated'] = dict(backward_error_normwise=nw, backward_error_componentwise=cw,
                                       backward_error_equilibrated=eta(xs),
                                       G=-float(xs[iV]) / DELTA_L1, r0=-float(xs[iV]) / DELTA_L1 / G - 1)
    gns = newton_step_G()
    out['G_NS_umfpack'] = dict(G=gns, r0=gns / G - 1)
    rec['result'] = out
    runlog.write_json(os.path.join(a.out, f'l1diag_L{a.L:g}_lev{a.level}_{a.mode}_VG{a.vg}.json'), rec)
    for s in steps:
        print('IR', s)
    print('EQ', out['superlu_equilibrated']); print('NS', out['G_NS_umfpack'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
