"""T1 version 2 (C-4 Task B): zero-drain-bias identity test with an equilibrated matrix-export
solver and an L2-native error estimate. PRE-REGISTERED: the definitions below were committed
before any T1-v2 result was produced (see phase3/specs/C4_T1v2_Preregistration_2026-10-30.md).
T1 v1 (phase3/c3/t1.py) and its results remain the historical record and are not re-evaluated.

Setting (unchanged from v1): V_D = 0 equilibrium converged to tau_eq = 1e-12 V (R05 D-3 rule);
G_GCA,h = (q mu_n / L) sum_i n0_i w_i from the independent 1COL solve on identical depth lines;
r0_X = G_X / G_GCA,h - 1. One linear system per state:
    A x = rhs,   A = exported static Newton matrix (bitwise the factored Jacobian, C-1(b) B1),
    rhs = (b(V1 = delta) - b(V1 = 0)) / delta,  delta = 1 mV  (one nonzero: the V1 circuit row),
    G(x) = c^T x with c = -e_{V1.I}.
AC (L2):  (A + j w A_d) x = rhs  (A_d = exported dynamic matrix; sign verified in development),
    G_L2 = Re G(x_AC(f2)), f1 = 1 Hz, f2 = 0.1 Hz.

Solutions computed per state:
    x_U   unscaled SciPy SuperLU of A (= the v1 L1 path; comparison only)
    x_E   equilibrated SuperLU (= L1 v2): A_s = D_r A D_c, D_r = 1/max_j |a_ij| (rows), then
          D_c = 1/max_i |(D_r A)_ij| (columns);  y = A_s^{-1} (D_r rhs),  x_E = D_c y.
          Consistent transformation of matrix, unknowns (x = D_c y) and right-hand side.
    x_M   UMFPACK 5.1 (DEVSIM's own solver, via devsim.umfpack.umfshim) on the exported (A, rhs)
    x_NS  DEVSIM's own first undamped Newton step at V1 = delta, divided by delta (as v1 G_NS)
    x_AC  DEVSIM small-signal AC solution vector at f2 (and f1 for E_freq)
Adjoint (output sensitivity):  A^T lam = c.
    lam_M via UMFPACK (transpose solve) -- used for L2 and for diagnostics;
    lam_E via the equilibrated SuperLU factorization (lam_E = D_r A_s^{-T} (D_c c)) -- used for L1v2.
Metrics for a solution x of system S (S = A, or A + j w A_d for AC), residual r = S x - rhs:
    eta(x)    normwise backward error of the row/column-equilibrated system (C-1(b) B3 metric):
              ||D_r r||_inf / (||A_s||_inf ||x/D_c||_inf + ||D_r rhs||_inf)
    E_dwr(x)  |lam^T r| / G_GCA,h   (dual-weighted residual: first-order estimate of G(x) - G*)
    s_i       = (|A||x| + |rhs|)_i : componentwise rounding scale of row i of the residual
    E_round(x) u ||lam o s||_2 / |G(x)|, u = 2^-53: probabilistic (root-sum-square) estimate of the
              rounding uncertainty of the residual, hence of E_dwr itself (Higham & Mary,
              SIAM J. Sci. Comput. 41(5), 2019: independent rounding errors grow like sqrt(n) u)
    E_cond(x) kappa_G u, kappa_G = |lam|^T s / |G(x)| (componentwise output condition number):
              worst-case first-order counterpart of E_round (all row errors aligned). REPORTED
              and used only for the conservative sensitivity verdict, not for the primary one.
    E_sec(x)  max_k |G_sec,k(x) - G(x)| / G_GCA,h, linearized SG current through mid-cell lines
              near L/4, L/2, 3L/4 of the channel, from the node values of the same x
E budgets (sums of ESTIMATES, never bounds):
    E_L1v2 = E_eq,L1v2 + E_sec(x_E) + E_dwr(x_E; lam_E) + E_round(x_E; lam_E)
    E_L2v2 = E_freq + E_eq,L2 + E_sec(x_AC) + E_dwr(x_AC; lam_M) + E_round(x_AC; lam_M)
    Conservative sensitivity budgets E_wc,X: the same with E_cond in place of E_round.
    E_eq,X = |r0_X(tau_eq) - r0_X(tau_eq/100)|;  E_freq = |ReY(f1) - ReY(f2)| / (0.99 G_GCA,h)
    No cross-path term: E_L2v2 does not use any L1 solution; E_L1v2 does not use x_AC or UMFPACK.
Acceptance (R04 8.1, UNCHANGED thresholds), applied to L2v2 (primary) and L1v2 separately:
    PASS if |r0| <= 3E and E <= 1e-6; NON-RESOLVING (class I) if E > 1e-6; FAIL otherwise.
    The same rule with E_wc gives the conservative sensitivity verdict (reported, not primary).
Diagnostics (reported, not part of the verdict):
    consistency |r0_L1v2 - r0_L2v2| <= 3 (E_L1v2 + E_L2v2);
    B3v2: eta(x_E) <= 1e-13;  B5v2: E_sec + E_dwr + E_round of x_E <= 1e-8;
    DWR validation: predicted G(x_U) - G(x_M) = lam_M^T r(x_U) vs actual;
    decoupling (R-I(a)): |G_E(s_couple=1) - G_E(s_couple=0)| / G_GCA,h <= E_L1v2;
    x_U reproduces the v1 G_L1 (same factorization, same rhs scaling as v1).
"""
from __future__ import annotations

import argparse
import array
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
from devsim.umfpack import umfshim as UMF

from phase3 import m0, runlog
from phase3.c1 import deck as D, measure as M
from phase3.c1.c1b_jacobian import export, snapshot, restore, to_vector
from phase3.c1.corner_study import eqnums, continuity, gauss_closure, oxide_lateral_rows
from phase3.c3.t1 import equilibrate, g_gca, onecol_sol, set_lateral_data
from phase3.c4 import deck4 as D4
from dataclasses import replace

TAU_EQ = 1e-12
DELTA = 1e-3
F1, F2 = 1.0, 0.1
U_ROUND = 2.0 ** -53
E_REQ = 1e-6
VERSION = 't1v2-2026-10-30'


# ----------------------------------------------------------------------------- linear algebra
def systems():
    A, b0, Ad, *_ = export()
    ds.circuit_alter(name='V1', value=DELTA)
    _, bd, *_ = export()
    ds.circuit_alter(name='V1', value=0.)
    return A.tocsc(), Ad.tocsc(), b0, bd


def equilibrate_matrix(A):
    dr = 1.0 / np.asarray(abs(A).max(axis=1).todense()).ravel()
    Ar = sp.diags(dr) @ A
    dc = 1.0 / np.asarray(abs(Ar).max(axis=0).todense()).ravel()
    As = (Ar @ sp.diags(dc)).tocsc()
    return As, dr, dc


class EqSolver:
    """Equilibrated SuperLU: A_s = D_r A D_c."""

    def __init__(self, A):
        self.As, self.dr, self.dc = equilibrate_matrix(A)
        self.lu = spla.splu(self.As)
        self.nAs = float(abs(self.As).sum(axis=1).max())

    def solve(self, b):
        return self.dc * self.lu.solve(self.dr * b)

    def solve_T(self, c):
        return self.dr * self.lu.solve(self.dc * c, trans='T')

    def eta(self, S, x, rhs):
        r = S @ x - rhs
        xe = x / self.dc; be = self.dr * rhs
        return float(np.max(np.abs(self.dr * r)) / (self.nAs * np.max(np.abs(xe)) + np.max(np.abs(be))))


def umf_solve(A, b, transpose=False):
    so = UMF.dsobject(A.shape[0], transpose)
    so.factor(Ap=array.array('i', A.indptr.astype(np.int32)), Ai=array.array('i', A.indices.astype(np.int32)),
              Ax=array.array('d', A.data), is_complex=False, is_same_symbolic=False)
    so.solve(b=array.array('d', b))
    return np.array(so.x)


def metrics(S, x, rhs, lam, absA, iV, G, eqs, geo, eq):
    r = S @ x - rhs
    Gx = -float(np.real(x[iV]))
    dwr = complex(lam @ r)
    sc = absA @ np.abs(x) + np.abs(rhs)
    kappa = float(np.abs(lam) @ sc) / abs(Gx)
    e_round = U_ROUND * float(np.linalg.norm(np.abs(lam) * sc)) / abs(Gx)
    secs = section_currents(geo, node_vectors(eq, x))
    gsec = [-float(np.real(s)) for s in secs]
    return dict(G=Gx, r0=Gx / G - 1, eta=eqs.eta(S, x, rhs), dwr_rel=float(np.real(dwr)) / G,
                dwr_imag_rel=float(np.imag(dwr)) / G, E_dwr=abs(float(np.real(dwr))) / G,
                kappa_G=kappa, E_cond=kappa * U_ROUND, E_round=e_round, G_sections=gsec,
                E_sec=max(abs(g - Gx) for g in gsec) / G, max_abs_residual=float(np.max(np.abs(r))))


# ----------------------------------------------------------------------------- state vectors
def node_vectors(eq, s):
    out = {}
    for v in ('Potential', 'Electrons'):
        rows = eq[(D.SI, v)]
        out[v] = np.where(rows >= 0, s[np.maximum(rows, 0)], 0.0)
    return out


def section_currents(geo, dv):
    e = M.edges(); out = []
    xs = geo.xs_ch; nch = len(xs) - 1; dx = xs[1] - xs[0]
    for kc in (round(nch / 4), nch // 2, round(3 * nch / 4) - 1):
        xc = float(xs[kc] + dx / 2)
        idx, sgn = M.crossing(xc)
        tot = 0j
        for v in ('Potential', 'Electrons'):
            for k, nn in (('n0', e['n0']), ('n1', e['n1'])):
                d = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name=f'ElectronCurrent:{v}@{k}'))
                tot += np.sum(sgn * d[idx] * dv[v][nn[idx]] * e['ec'][idx])
        out.append(tot)
    return out


def newton_step_vector(eq, n):
    X0 = snapshot()
    D.redeclare_updates('default', 'default')
    ds.circuit_alter(name='V1', value=DELTA)
    D.solve(abs_err=1e300, rel_err=1e300, it=1)
    X1 = snapshot()
    restore(X0); ds.circuit_alter(name='V1', value=0.)
    D.redeclare_updates('log_damp', 'positive')
    dX = {k: (X1[k] - X0[k]) / DELTA for k in X0}
    return to_vector(dX, eq, n)


def ac_vector(eq, n, f):
    ds.solve(type='ac', frequency=f, solver_type='direct')
    sr, si = {}, {}
    for r, v in ((D.SI, 'Potential'), (D.SI, 'Electrons'), (D.OX, 'Potential')):
        sr[(r, v)] = np.array(ds.get_node_model_values(device=D.DEV, region=r, name=f'{v}_Real'))
        si[(r, v)] = np.array(ds.get_node_model_values(device=D.DEV, region=r, name=f'{v}_Imag'))
    for nd in ds.get_circuit_node_list():
        sr[('circuit', nd)] = np.array([ds.get_circuit_node_value(solution='ssac_real', node=nd)])
        si[('circuit', nd)] = np.array([ds.get_circuit_node_value(solution='ssac_imag', node=nd)])
    return to_vector(sr, eq, n) + 1j * to_vector(si, eq, n)


# ----------------------------------------------------------------------------- one evaluation
def evaluate_light(G):
    """Tightened-equilibrium repeat for E_eq: only the two quantities E_eq needs
    (G of the equilibrated SuperLU solution and Re of the AC solution at f2)."""
    eq = eqnums()
    A, Ad, b0, bd = systems()
    n = A.shape[0]; iV = int(eq[('circuit', 'V1.I')][0])
    x_E = EqSolver(A).solve((bd - b0) / DELTA)
    x_AC = ac_vector(eq, n, F2)
    return dict(E=dict(r0=-float(x_E[iV]) / G - 1), AC_f2=dict(r0=-float(np.real(x_AC[iV])) / G - 1))


def evaluate(geo, G, compare_unscaled):
    """All solutions and metrics at the current (converged) state."""
    eq = eqnums()
    A, Ad, b0, bd = systems()
    n = A.shape[0]; iV = int(eq[('circuit', 'V1.I')][0])
    rhs = (bd - b0) / DELTA
    absA = abs(A)
    eqs = EqSolver(A)
    c = np.zeros(n); c[iV] = -1.0
    lam_E = eqs.solve_T(c)
    lam_M = umf_solve(A, c, transpose=True)
    out = dict(n_unknowns=int(n), rhs_nnz=int(np.count_nonzero(bd - b0)),
               adjoint=dict(res_lam_E=float(np.max(np.abs(A.T @ lam_E - c))),
                            res_lam_M=float(np.max(np.abs(A.T @ lam_M - c))),
                            lam_rel_diff=float(np.max(np.abs(lam_E - lam_M)) / np.max(np.abs(lam_M)))))
    x_E = eqs.solve(rhs)
    out['E'] = metrics(A, x_E, rhs, lam_E, absA, iV, G, eqs, geo, eq)
    out['E']['dwr_with_lam_M_rel'] = float(lam_M @ (A @ x_E - rhs)) / G
    x_M = umf_solve(A, rhs)
    out['M'] = metrics(A, x_M, rhs, lam_M, absA, iV, G, eqs, geo, eq)
    if compare_unscaled:
        lu = spla.splu(A)
        x_U = lu.solve(bd - b0) / DELTA                       # exactly the v1 L1 expression
        out['U'] = metrics(A, x_U, rhs, lam_M, absA, iV, G, eqs, geo, eq)
        out['U']['dwr_with_lam_E_rel'] = float(lam_E @ (A @ x_U - rhs)) / G
        out['U']['actual_minus_M_rel'] = (out['U']['G'] - out['M']['G']) / G
    x_NS = newton_step_vector(eq, n)
    out['NS'] = metrics(A, x_NS, rhs, lam_M, absA, iV, G, eqs, geo, eq)
    ac = {}
    for f in (F1, F2):
        x_AC = ac_vector(eq, n, f)
        S = (A + 1j * 2 * math.pi * f * Ad).tocsc()
        ac[f] = metrics(S, x_AC, rhs, lam_M, absA, iV, G, eqs, geo, eq)
        ac[f]['ImY_rel'] = float(np.imag(-x_AC[iV])) / G
    out['AC_f1'], out['AC_f2'] = ac[F1], ac[F2]
    return out


def bottom_rows(geo, eq):
    """Equations at the bottom-boundary nodes of the lateral/strip contacts and next to them."""
    A, *_ = export(); Ar = A.tocsr()
    xs, ys = map(np.array, D.coords(D.SI))
    psi = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Potential'))
    n = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Electrons'))
    ymax = float(ys.max()); out = []
    picks = sorted({0, geo.ic0, min(geo.ic0 + 1, len(geo.xs) - 1), geo.ic1 - 1, geo.ic1, len(geo.xs) - 1})
    for ic in picks:
        nd = int(np.argmin((xs - geo.xs[ic]) ** 2 + (ys - ymax) ** 2))
        row = {}
        for v in ('Potential', 'Electrons'):
            rr = int(eq[(D.SI, v)][nd]); row[f'{v}_row_nnz'] = int(Ar.indptr[rr + 1] - Ar.indptr[rr])
        out.append(dict(x_cm=float(xs[nd]), column=ic, psi_minus_psib0=float(psi[nd] - D.PSI_B0),
                        phin_V=float(psi[nd] - D.VT * math.log(n[nd] / D.NI)), **row))
    return out


def t1_point(geo, vg, model, compare_unscaled):
    mode = D4.mode_of(model)
    r = dict(vg=vg, model=model, mode=mode)
    ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
    sol0, r['equilibrium'] = equilibrate(geo, vg, mode, TAU_EQ)
    G, ns = g_gca(geo, sol0)
    r.update(G_GCA_h_S_per_cm=G, Ns_FV_cm2=ns, G_2D_midcolumn_S_per_cm=M.gca_conductance(geo.L_cm)['G_S_per_cm'])
    ev = evaluate(geo, G, compare_unscaled); r['eval'] = ev
    # tightened equilibrium for E_eq (both paths, recomputed from scratch)
    rt = D.newton(TAU_EQ / 100, extra=0)
    Gt, _ = g_gca(geo, onecol_sol(geo, vg, tol=TAU_EQ / 100))
    ev_t = evaluate_light(Gt)
    r['tight'] = dict(newton_iterations=rt['iterations'], G_GCA_h=Gt,
                      r0_E=ev_t['E']['r0'], r0_AC_f2=ev_t['AC_f2']['r0'])
    E_eq_L1 = abs(ev_t['E']['r0'] - ev['E']['r0']); E_eq_L2 = abs(ev_t['AC_f2']['r0'] - ev['AC_f2']['r0'])
    E_freq = abs(ev['AC_f1']['G'] - ev['AC_f2']['G']) / (0.99 * G)
    e, a = ev['E'], ev['AC_f2']
    L1 = dict(r0=e['r0'], E_eq=E_eq_L1, E_sec=e['E_sec'], E_dwr=e['E_dwr'], E_round=e['E_round'], E_cond=e['E_cond'])
    L1['E'] = E_eq_L1 + e['E_sec'] + e['E_dwr'] + e['E_round']
    L1['E_wc'] = E_eq_L1 + e['E_sec'] + e['E_dwr'] + e['E_cond']
    L2 = dict(r0=a['r0'], E_freq=E_freq, E_eq=E_eq_L2, E_sec=a['E_sec'], E_dwr=a['E_dwr'], E_round=a['E_round'], E_cond=a['E_cond'])
    L2['E'] = E_freq + E_eq_L2 + a['E_sec'] + a['E_dwr'] + a['E_round']
    L2['E_wc'] = E_freq + E_eq_L2 + a['E_sec'] + a['E_dwr'] + a['E_cond']

    def verdict(r0, E):
        return 'NON-RESOLVING (class I)' if E > E_REQ else 'PASS' if abs(r0) <= 3 * E else 'FAIL'
    for X in (L1, L2):
        X['verdict'] = verdict(X['r0'], X['E'])
        X['verdict_wc'] = verdict(X['r0'], X['E_wc'])
    r['L1v2'], r['L2v2'] = L1, L2
    r['diagnostics'] = dict(
        consistency=dict(diff=abs(L1['r0'] - L2['r0']), limit=3 * (L1['E'] + L2['E']),
                         ok=abs(L1['r0'] - L2['r0']) <= 3 * (L1['E'] + L2['E'])),
        B3v2=dict(eta_E=e['eta'], ok=e['eta'] <= 1e-13),
        B5v2=dict(lin=e['E_sec'] + e['E_dwr'] + e['E_round'], ok=e['E_sec'] + e['E_dwr'] + e['E_round'] <= 1e-8),
        NS_vs_M_rel=(ev['NS']['G'] - ev['M']['G']) / G, AC_vs_M_rel=(a['G'] - ev['M']['G']) / G)
    if compare_unscaled:
        u = ev['U']
        r['diagnostics']['dwr_validation_U'] = dict(predicted_lam_M=u['dwr_rel'], predicted_lam_E=u['dwr_with_lam_E_rel'],
                                                     actual_U_minus_M=u['actual_minus_M_rel'])
    if mode == 'a':
        ds.set_parameter(device=D.DEV, name='s_couple', value=1.)
        D.converge()
        eqn = eqnums(); A, _, b0, bd = systems(); iV = int(eqn[('circuit', 'V1.I')][0])
        gc = -float(EqSolver(A).solve((bd - b0) / DELTA)[iV])
        r['decoupling'] = dict(G_E_s1=gc, rel_change=(gc - e['G']) / G, ok=abs(gc - e['G']) / G <= L1['E'])
        ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
        D.converge()
    eq = eqnums(); A, *_ = export()
    r['continuity'] = continuity()
    r['gauss'] = gauss_closure(eq)
    r['oxide_lateral_single_entry_rows'] = oxide_lateral_rows(eq, A)['single_entry_rows']
    r['uniformity'] = [M.lateral_uniformity('Potential'), M.lateral_uniformity('Electrons')]
    r['bottom_rows'] = bottom_rows(geo, eq)
    r['contacts'] = {c: ds.get_contact_equation_list(device=D.DEV, contact=c) for c in ds.get_contact_list(device=D.DEV)}
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True, choices=['a', 'b_c3', 'b', 'c'])
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--lres', type=float, default=0.0, help='l_res / W_dep (model c only)')
    ap.add_argument('--vgs', default='0.8,1.25,2.0')
    ap.add_argument('--compare-unscaled', type=int, default=None)
    ap.add_argument('--out', required=True)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    cu = a.compare_unscaled if a.compare_unscaled is not None else int(a.model in ('a', 'b_c3'))
    geo = D4.Geo4(a.model, a.L, a.level, a.lres)
    rec = runlog.header(dict(version=VERSION, model=a.model, L_um=a.L, level=a.level, lres=a.lres, vgs=a.vgs,
                             compare_unscaled=cu, tau_eq=TAU_EQ, delta=DELTA, f1=F1, f2=F2, E_req=E_REQ,
                             pm=D.pm_checklist()))
    t0 = time.time()
    rec['mesh'] = D4.build(geo)
    P = replace(m0.M0Params.f0(), L=a.L * 1e-6)
    rec['points'] = {}
    for vg in [float(v) for v in a.vgs.split(',')]:
        t = time.time()
        p = t1_point(geo, vg, a.model, bool(cu))
        p['M0_G0_S_per_cm'] = m0.conductance0(vg, P) * 1e4
        p['wall_s'] = time.time() - t
        rec['points'][f'{a.model}:VG={vg}'] = p
        print('T1v2', a.model, a.lres, a.L, a.level, vg, json.dumps(dict(L1v2=p['L1v2'], L2v2=p['L2v2'])), flush=True)
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    tag = f'{a.model}' + (f'_lres{a.lres:g}' if a.model == 'c' else '')
    runlog.write_json(os.path.join(out, f't1v2_{tag}_L{a.L:g}_lev{a.level}.json'), rec)
    return 0


if __name__ == '__main__':
    sys.exit(main())
