"""C-1(b) numerical acceptance study (Task C, 2026-10-30).

Same deck and states as C-1(b) (legacy built-in-mesher deck, V_G = 1.25 V, V_D = 0
converged, then V1 stepped to delta). The ORIGINAL observation is retained: the
external SuperLU update differed from DEVSIM's UMFPACK update by 1.4-1.9e-10.

Measured here, per state:
  1. bitwise reproduction of the stored C-1(b) system (repro_f4ca0bf npz);
  2. residual norms of each solution (raw, per block);
  3. backward errors: normwise raw, normwise after row/column equilibration,
     componentwise (Oettli-Prager) over all rows and over rows that are not
     structurally decoupled (single-entry rows with b = 0 excluded);
  4. condition estimates (Hager-Higham 1-norm estimator = LOWER-bound estimate):
     raw, row-equilibrated, row+column-equilibrated; Skeel componentwise condition
     per unknown block cond_B(A,x) = || (|A^-1| |A| |x|)_B ||_inf / ||x_B||_inf;
  5. a reference solution by mixed-precision iterative refinement (residuals in
     80-bit long double), and the forward error of each double-precision solver
     against it, per block; double-precision refinement of the SuperLU solution;
  6. directional-derivative checks of the Jacobian: central differences of the
     exported residual F = -rhs along random, variable-scaled directions over a
     ladder of step sizes (truncation/round-off V-shape);
  7. damped vs undamped Newton steps: DEVSIM 'log_damp'/'positive' vs 'default',
     with the damping reproduced from the undamped update using the formulas of
     DEVSIM Equation.cc (LogSolutionUpdate, PositiveSolutionUpdate).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import devsim as ds

from phase3 import runlog
from phase3.c1 import deck as D
from phase3.c1.c1b_jacobian import REC, recording_callback, export, snapshot, restore, to_vector

REPRO = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                     'results/phase3/C1/repro_f4ca0bf/c1b')


def eqnums():
    stock = ds.get_parameter(name='solver_callback')
    REC.clear(); X0 = snapshot()
    ds.set_parameter(name='direct_solver', value='custom')
    ds.set_parameter(name='solver_callback', value=recording_callback)
    try:
        D.solve(abs_err=1e300, rel_err=1e300, it=1)
    finally:
        ds.set_parameter(name='solver_callback', value=stock)
    restore(X0)
    return REC['eqnums']


def blocks(eq):
    return {f'{k[0]}:{k[1]}': v[v >= 0] for k, v in eq.items() if len(v[v >= 0])}


def builtin_step():
    X0 = snapshot(); D.solve(abs_err=1e300, rel_err=1e300, it=1); X1 = snapshot(); restore(X0)
    return X0, X1


def csr_matvec_ld(A, x):
    """y = A x with products and sums in numpy long double (80-bit on x86-64)."""
    A = A.tocsr()
    data = A.data.astype(np.longdouble); xl = np.asarray(x, dtype=np.longdouble)
    prod = data * xl[A.indices]
    y = np.zeros(A.shape[0], dtype=np.longdouble)
    nz = np.diff(A.indptr) > 0
    y[nz] = np.add.reduceat(prod, A.indptr[:-1][nz])
    return y


def analyse_system(A, b, xs, eq):
    n = A.shape[0]; out = {}
    absA = abs(A); B = blocks(eq)
    rowmax = np.asarray(absA.max(axis=1).todense()).ravel()
    Dr = sp.diags(1 / rowmax); As = (Dr @ A).tocsc()
    colmax = np.asarray(abs(As).max(axis=0).todense()).ravel()
    Dc = sp.diags(1 / colmax); Ae = (As @ Dc).tocsc(); be = Dr @ b
    single = np.diff(A.tocsr().indptr) == 1
    decoupled = single & (b == 0)
    out['structure'] = dict(n=n, nnz=int(A.nnz), single_entry_rows=int(single.sum()),
                            decoupled_rows=int(decoupled.sum()))
    nA_inf = float(abs(A).sum(axis=1).max()); nAe_inf = float(abs(Ae).sum(axis=1).max())
    for name, x in xs.items():
        r = A @ x - b
        den = absA @ np.abs(x) + np.abs(b)
        m = den > 0; om_all = float(np.max(np.abs(r[m]) / den[m]))
        m2 = m & ~decoupled; om_nd = float(np.max(np.abs(r[m2]) / den[m2]))
        xe = x / Dc.diagonal(); re = Ae @ xe - be
        wi = np.where(m2)[0][int(np.argmax(np.abs(r[m2]) / den[m2]))]
        wblk = next((k for k, v in B.items() if wi in set(v.tolist())), None)
        out[name] = dict(
            omega_dominant_row=dict(row=int(wi), block=wblk, row_nnz=int(np.diff(A.tocsr().indptr)[wi]),
                                    b=float(b[wi]), x=float(x[wi]), den=float(den[wi]),
                                    den_over_max_den=float(den[wi] / den.max())),
            residual_inf=float(np.max(np.abs(r))), residual_2=float(np.linalg.norm(r)),
            residual_inf_per_block={k: float(np.max(np.abs(r[v]))) for k, v in B.items()},
            eta_normwise_raw=float(np.max(np.abs(r)) / (nA_inf * np.max(np.abs(x)) + np.max(np.abs(b)))),
            eta_normwise_equilibrated=float(np.max(np.abs(re)) / (nAe_inf * np.max(np.abs(xe)) + np.max(np.abs(be)))),
            omega_componentwise_all_rows=om_all, omega_componentwise_excl_decoupled=om_nd)

    def cond1(M):
        lu = spla.splu(M.tocsc())
        Minv = spla.LinearOperator(M.shape, matvec=lu.solve, rmatvec=lambda v: lu.solve(v, trans='T'))
        return float(spla.onenormest(M) * spla.onenormest(Minv))
    out['cond1_estimates_lower_bound'] = dict(raw=cond1(A), row_equilibrated=cond1(As),
                                              row_col_equilibrated=cond1(Ae))
    # Skeel per block: || P_B |A^-1| g ||_inf / ||x_B||_inf,  g = |A||x|
    lu = spla.splu(A.tocsc())
    x = xs['builtin_umfpack']; g = absA @ np.abs(x)
    sk = {}
    for k, idx in list(B.items()) + [('ALL', np.arange(n))]:
        mask = np.zeros(n); mask[idx] = 1.0
        xb = np.max(np.abs(x[idx]))
        # C = mask * (A^-1 diag(g)) (n x n);  ||C||_inf = ||C^T||_1, C^T = diag(g) A^-T mask
        op = spla.LinearOperator((n, n), matvec=lambda v, m=mask: g * lu.solve(m * np.ravel(v), trans='T'),
                                 rmatvec=lambda v, m=mask: m * lu.solve(g * np.ravel(v)))
        try:
            sk[k] = float(spla.onenormest(op) / xb) if xb > 0 else None
        except Exception as exc:
            sk[k] = f'failed {exc!r}'
    out['skeel_cond_per_block_estimate'] = sk
    return out, lu


def refine_reference(A, b, lu, x0, iters=8):
    x = np.asarray(x0, dtype=np.longdouble); hist = []
    bl = np.asarray(b, dtype=np.longdouble)
    for _ in range(iters):
        r = bl - csr_matvec_ld(A, x)
        dx = lu.solve(np.asarray(r, dtype=float))
        x = x + np.asarray(dx, dtype=np.longdouble)
        hist.append(float(np.max(np.abs(dx)) / np.max(np.abs(np.asarray(x, float)))))
    return np.asarray(x, dtype=float), hist


def refine_double(A, b, lu, x0, iters=3):
    x = x0.copy(); hist = []
    for _ in range(iters):
        dx = lu.solve(b - A @ x); x = x + dx
        hist.append(float(np.max(np.abs(dx)) / np.max(np.abs(x))))
    return x, hist


def blockrel(a, ref, B):
    return {k: float(np.max(np.abs(a[v] - ref[v])) / np.max(np.abs(ref[v]))) if np.max(np.abs(ref[v])) > 0 else None
            for k, v in B.items()}


def directional(eq, X0, A, b, n_dirs=4, hs=(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)):
    """F = -rhs. Directions v scaled per unknown: psi 1 V, n |n|, circuit 1 (V or A/cm
    scaled by |value| + 1e-6). Step t*h with t = 1; compares (F(x+hv)-F(x-hv))/2h to A v."""
    rng = np.random.default_rng(2026)
    n = A.shape[0]; res = []
    keys = list(X0.keys())
    for d in range(n_dirs):
        v = {}
        for k in keys:
            z = rng.standard_normal(len(X0[k]))
            if k[1] == 'Electrons':
                v[k] = z * np.abs(X0[k])
            elif k[0] == 'circuit':
                v[k] = z * (np.abs(X0[k]) + 1e-6)
            else:
                v[k] = z * 1.0
        vv = to_vector(v, eq, n); Av = A @ vv
        row = []
        for h in hs:
            F = []
            for s in (+1, -1):
                restore({k: X0[k] + s * h * v[k] for k in keys})
                F.append(-export()[1])
            fd = (F[0] - F[1]) / (2 * h)
            row.append(dict(h=h, rel_err=float(np.max(np.abs(fd - Av)) / np.max(np.abs(Av)))))
        restore(X0)
        res.append(dict(direction=d, ladder=row, best=min(row, key=lambda q: q['rel_err'])))
    return res


def damping_study(eq, delta, outdir):
    n = None
    ds.circuit_alter(name='V1', value=delta)
    X0 = snapshot()
    D.redeclare_updates('default', 'default'); _, Xu = builtin_step()
    D.redeclare_updates('log_damp', 'positive'); _, Xd = builtin_step()
    D.redeclare_updates('default', 'default')
    out = dict(delta_V=delta, blocks={})
    for k in X0:
        du = Xu[k] - X0[k]; dd = Xd[k] - X0[k]
        if k[1] == 'Potential':
            pred = np.where(np.abs(du) > 0.0259, np.sign(du) * 0.0259 * np.log1p(np.abs(du) / 0.0259), du)
        elif k[1] == 'Electrons':
            nv = X0[k] + du
            pred = np.where(nv <= 0, 0.001 * X0[k] - X0[k], du)
        else:
            pred = du
        den = np.max(np.abs(du)) if np.max(np.abs(du)) > 0 else 1.0
        out['blocks'][f'{k[0]}:{k[1]}'] = dict(
            max_abs_undamped=float(np.max(np.abs(du))),
            n_components_damped=int(np.sum(np.abs(dd - du) > 1e-15 * max(1.0, den))),
            damped_vs_undamped_rel=float(np.max(np.abs(dd - du)) / den),
            damped_vs_formula_rel=float(np.max(np.abs(dd - pred)) / den))
    restore(X0)
    ds.circuit_alter(name='V1', value=0.0)
    D.converge()
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=False)
    rec = runlog.header(dict(test='C-1(b) acceptance study', pm=D.pm_checklist()))
    t0 = time.time()
    rec['mesh'] = D.build_mesh(); D.build_physics()
    rec['equilibrium'] = D.equilibrium(1.25)
    D.redeclare_updates('default', 'default')
    rec['undamped_reconverge'] = D.converge()
    rec['states'] = {}
    for tag, delta in (('dV1mV', 1e-3), ('dV50mV', 5e-2)):
        s = {}
        ds.circuit_alter(name='V1', value=delta)
        X0 = snapshot(); eq = eqnums(); n_unk = sum(len(v) for v in eq.values())
        A, b, *_ = export()
        st = np.load(os.path.join(REPRO, f'{tag}_system.npz'))
        As = sp.csc_matrix((st['A_data'], st['A_indices'], st['A_indptr']), shape=A.shape)
        s['reproduces_stored_system'] = dict(A_max_abs_diff=float(abs(A - As).max()) if (A - As).nnz else 0.0,
                                             rhs_max_abs_diff=float(np.max(np.abs(b - st['rhs']))))
        _, X1 = builtin_step()
        x_bi = to_vector({k: X1[k] - X0[k] for k in X0}, eq, A.shape[0])
        x_ext = spla.splu(A.tocsc()).solve(b)
        s['original_observation'] = dict(stored_dx_builtin_vs_now=float(np.max(np.abs(x_bi - st['dx_builtin']))),
                                         cross_solver_rel_per_block=blockrel(x_ext, x_bi, blocks(eq)))
        ana, lu = analyse_system(A, b, dict(builtin_umfpack=x_bi, scipy_superlu=x_ext), eq)
        s['analysis'] = ana
        x_ref, hist = refine_reference(A, b, lu, x_ext)
        x_dr, hist_d = refine_double(A, b, lu, x_ext)
        B = blocks(eq)
        s['reference_refinement'] = dict(longdouble_eps=float(np.finfo(np.longdouble).eps),
                                         update_history=hist,
                                         fwd_err_builtin=blockrel(x_bi, x_ref, B),
                                         fwd_err_superlu=blockrel(x_ext, x_ref, B),
                                         fwd_err_superlu_double_refined=blockrel(x_dr, x_ref, B),
                                         double_refinement_history=hist_d,
                                         refined_vs_builtin=blockrel(x_dr, x_bi, B))
        s['directional_derivative'] = directional(eq, X0, A, b)
        restore(X0)
        ds.circuit_alter(name='V1', value=0.0); D.converge()
        rec['states'][tag] = s
    rec['damping'] = [damping_study(eqnums(), d, out) for d in (5e-2, 2e-1)]
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(out, 'c1b_acceptance.json'), rec)
    print('C1B_ACC_DONE')
    return 0


if __name__ == '__main__':
    sys.exit(main())
