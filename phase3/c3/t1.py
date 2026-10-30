"""Formal zero-drain-bias consistency test T1 for R-I(a) and R-I(b).

Specification: Research_05 v2.1 C-3/C-4, Research_04 v2.1 5.2 and 8.1, Research_07.
Deck: explicit tensor mesh with the C-2 triple-point continuity tie (deck.build_mesh_explicit
corner='tie'); D-1 depth/oxide maps; lateral map uniform N_y,0 = 50; levels k = 0..3.

PRE-REGISTERED DEFINITIONS (fixed before any result was produced):
  G_GCA,h  = (q mu_n / L) * sum_i n0_i w_i, n0 from the independent 1COL FV solve on the
             identical depth lines (onecol.py), w_i = 1-D dual widths (D4).  [S/cm]
  G_L1     = -dI_V1/dV1 from A^-1 (b(V1=delta) - b(V1=0)), exported static Newton matrix,
             SciPy SuperLU, delta = 1 mV (b is exactly linear in delta).
  G_NS     = same quantity from DEVSIM's own first undamped Newton step (UMFPACK).
  G_L2(f)  = Re Y(f) = -Re I_V1(f) for a 1 V AC source, DEVSIM small-signal AC.
  r0_X     = G_X / G_GCA,h - 1.
  E_eq     = |r0(tau_eq) - r0(tau_eq/100)|, tau_eq = 1e-12 V (2-D Newton and 1COL).
  E_freq   = |ReY(f1) - ReY(f2)| / (0.99 G_GCA,h), f1 = 1 Hz, f2 = 0.1 Hz (L2 only).
  E_lin,L1 = max_k |G_sec,k - G_L1| / G_GCA,h  +  |G_L1 - G_NS| / G_GCA,h
             (contact current vs linearized SG cross-section sums at mid-cell lines near
             L/4, L/2, 3L/4 of the same linear solution, plus cross-solver sensitivity).
  E_lin,L2 = max_k |Re I_sec,k(f2) - ReY(f2)| / G_GCA,h  +  |ReY(f2) - G_L1| / G_GCA,h.
  E_X      = E_freq + E_eq + E_lin,X   (sum of ESTIMATES; not a bound).
  Acceptance (R04 8.1, unchanged): |r0| <= 3E and E <= 1e-6, on every level and L.
     E > 1e-6 -> non-resolving (class I).  |r0| > 3E with E <= 1e-6 -> FAIL (N or P).
  Consistency: |r0_primary - r0_L3| <= 3 (E_primary + E_L3) when L3 was run.
  Decoupling (R-I(a) only): |G_L1(s_couple=1) - G_L1(s_couple=0)| / G_GCA,h <= E_L1.
  L3: delta_k = 0.25*2^k mV, k = 0..5; 2-D Newton tau = 1e-13 V (repeat at 1e-15 for e_i);
      R-I(a) drain psi data = 1COL(phi_n = delta); GCA side 4-point Gauss-Legendre over
      1COL slices (3-point for the quadrature estimate); Neville tableau, T and Lambda_m
      E_noise as R04 5.2; (k, m) with T monotone decreasing in m, minimising T + E_noise.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import scipy.sparse.linalg as spla

import devsim as ds

from phase3 import m0, meshmap as MM, runlog
from phase3.c1 import deck as D, measure as M, onecol
from phase3.c1.c1b_jacobian import export, snapshot, restore
from phase3.c1.corner_study import eqnums, continuity, gauss_closure, oxide_lateral_rows

TAU_EQ = 1e-12
DELTA_L1 = 1e-3
F1, F2 = 1.0, 0.1
L3_DELTAS = [0.25e-3 * 2 ** k for k in range(6)]
LAMBDA = {1: 3.0, 2: 5.0, 3: 6.428571428571429}       # R04 5.2 (ratio-2 nodes, to 0)


# ----------------------------------------------------------------------------- setup
class Geo:
    def __init__(self, L_um, level, n_lat0=50):
        self.xs = MM.lateral_map(L_um, n_lat0).nodes(level)
        self.y_si = MM.f0_depth_map().nodes(level)
        self.y_ox = MM.oxide_map(20.0, 5).nodes(level)
        self.L_cm = float(self.xs[-1]); self.L_um = L_um; self.level = level

    def quality(self):
        return dict(lateral=MM.lateral_map(self.L_um, 50).quality(self.level),
                    depth=MM.f0_depth_map().quality(self.level),
                    oxide=MM.oxide_map(20.0, 5).quality(self.level))


def onecol_sol(geo, vg, phin=0.0, tol=TAU_EQ):
    return onecol.solve_1col(geo.y_ox, geo.y_si, vg, tol=tol, phin=phin)


def g_gca(geo, sol):
    ns = float(np.sum(sol['n'] * sol['w_si']))
    return D.Q * D.MU_N * ns / geo.L_cm, ns


def _lut(sol):
    return {round(float(y), 15): float(p) for y, p in zip(sol['y'], sol['psi'])}


def set_lateral_data(geo, vg, sol0, sol_drain=None):
    """PsiLat = 1COL(phi_n = 0) everywhere; drain column (x = L) = 1COL(phi_n = V_D)."""
    xs, ys = map(np.array, D.coords(D.SI)); xo, yo = map(np.array, D.coords(D.OX))
    l0 = _lut(sol0); ld = _lut(sol_drain) if sol_drain is not None else l0
    L = geo.L_cm
    vs = [ld[round(float(y), 15)] if x == L else l0[round(float(y), 15)] for x, y in zip(xs, ys)]
    vo = [ld[round(float(y), 15)] if x == L else l0[round(float(y), 15)] for x, y in zip(xo, yo)]
    ds.set_node_values(device=D.DEV, region=D.SI, name='PsiLat', values=vs)
    ds.set_node_values(device=D.DEV, region=D.OX, name='PsiLat', values=vo)


def seed(sol0):
    l0 = _lut(sol0)
    for r in (D.SI, D.OX):
        _, ys = D.coords(r)
        ds.set_node_values(device=D.DEV, region=r, name='Potential', values=[l0[round(float(y), 15)] for y in ys])
    psi = ds.get_node_model_values(device=D.DEV, region=D.SI, name='Potential')
    ds.set_node_values(device=D.DEV, region=D.SI, name='Electrons', values=[D.NI * math.exp(p / D.VT) for p in psi])


def equilibrate(geo, vg, mode, tau=TAU_EQ):
    sol0 = onecol_sol(geo, vg)
    ds.set_parameter(device=D.DEV, name='gate_bias', value=vg)
    ds.circuit_alter(name='V1', value=0.)
    set_lateral_data(geo, vg, sol0); seed(sol0)
    D.set_lateral_dirichlet(mode == 'a')
    r = D.newton(tau, extra=0)
    if not r['converged']:
        raise D.NotConverged(f'equilibrate {mode} {vg}: {r["history"][-3:]}')
    return sol0, dict(onecol_iterations=sol0['iterations'], onecol_last_update=sol0['last_update_V'],
                      newton_iterations=r['iterations'], newton_history=r['history'])


# ----------------------------------------------------------------------------- linear response
def node_vectors(eq, s):
    out = {}
    for v in ('Potential', 'Electrons'):
        rows = eq[(D.SI, v)]
        out[v] = np.where(rows >= 0, s[np.maximum(rows, 0)], 0.0)
    return out


def section_currents(geo, dv):
    """Linearized SG electron current through mid-cell lines near L/4, L/2, 3L/4 for a
    node perturbation dv = {'Potential': array, 'Electrons': array} (may be complex)."""
    e = M.edges(); out = []
    nx = len(geo.xs); dx = geo.xs[1] - geo.xs[0]
    for kc in (round((nx - 1) / 4), (nx - 1) // 2, round(3 * (nx - 1) / 4) - 1):
        xc = float(geo.xs[kc] + dx / 2)
        idx, sgn = M.crossing(xc)
        tot = 0j
        for v in ('Potential', 'Electrons'):
            for k, nn in (('n0', e['n0']), ('n1', e['n1'])):
                d = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name=f'ElectronCurrent:{v}@{k}'))
                tot += np.sum(sgn * d[idx] * dv[v][nn[idx]] * e['ec'][idx])
        out.append(dict(x_cm=xc, I=tot))
    return out


def l1(geo, eq):
    A, b0, *_ = export()
    lu = spla.splu(A.tocsc())
    ds.circuit_alter(name='V1', value=DELTA_L1)
    _, bd, *_ = export()
    ds.circuit_alter(name='V1', value=0.)
    s = lu.solve(bd - b0) / DELTA_L1
    G = -float(s[int(eq[('circuit', 'V1.I')][0])])
    # section sums of the same linear solution; sign: section sum = -I_drain (C-1(c))
    secs = section_currents(geo, node_vectors(eq, s))
    gsec = [-float(x['I'].real) for x in secs]
    return G, gsec, int(A.shape[0])


def newton_step_G():
    X0 = snapshot()
    D.redeclare_updates('default', 'default')
    ds.circuit_alter(name='V1', value=DELTA_L1)
    D.solve(abs_err=1e300, rel_err=1e300, it=1)
    X1 = snapshot()
    restore(X0); ds.circuit_alter(name='V1', value=0.)
    D.redeclare_updates('log_damp', 'positive')
    return -(X1[('circuit', 'V1.I')][0] - X0[('circuit', 'V1.I')][0]) / DELTA_L1


def l2(geo, f):
    ds.solve(type='ac', frequency=f, solver_type='direct')
    Y = -(ds.get_circuit_node_value(solution='ssac_real', node='V1.I')
          + 1j * ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
    dv = {v: np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name=f'{v}_Real'))
          + 1j * np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name=f'{v}_Imag'))
          for v in ('Potential', 'Electrons')}
    secs = section_currents(geo, dv)
    return Y, [-x['I'] for x in secs]


# ----------------------------------------------------------------------------- L3
def gl_nodes(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def i_gca(geo, vg, delta, n):
    x, w = gl_nodes(n)
    phis = 0.5 * delta * (x + 1); tot = 0.0
    for p, ww in zip(phis, w):
        tot += ww * g_gca(geo, onecol_sol(geo, vg, phin=float(p)))[1]
    return D.Q * D.MU_N / geo.L_cm * 0.5 * delta * tot


def drain_current_at(geo, vg, mode, delta, sol0, tau):
    if mode == 'a':
        set_lateral_data(geo, vg, sol0, onecol_sol(geo, vg, phin=delta))
    ds.circuit_alter(name='V1', value=delta)
    r = D.newton(tau, extra=0, max_it=60)
    if not r['converged']:
        raise D.NotConverged(f'L3 {mode} {vg} {delta}')
    return -ds.get_circuit_node_value(solution='dcop', node='V1.I'), r['iterations']


def neville(ds_, rs, e):
    """R_{k,m}: degree-m polynomial through delta_k..delta_{k+m}, evaluated at 0."""
    n = len(ds_); R = {}
    for k in range(n):
        R[(k, 0)] = rs[k]
    for m in range(1, n):
        for k in range(n - m):
            dk, dkm = ds_[k], ds_[k + m]
            R[(k, m)] = (dkm * R[(k, m - 1)] - dk * R[(k + 1, m - 1)]) / (dkm - dk)
    best = None; table = []
    emax = max(e)
    for k in range(n):
        Ts = []
        for m in range(1, min(4, n - k)):
            T = abs(R[(k, m)] - R[(k, m - 1)])
            Ts.append(T)
            mono = all(Ts[i + 1] < Ts[i] for i in range(len(Ts) - 1))
            En = LAMBDA[m] * emax
            table.append(dict(k=k, m=m, R=R[(k, m)], T=T, E_noise=En, T_monotone=mono))
            if mono and m >= 2 and (best is None or T + En < best['T'] + best['E_noise']):
                best = dict(k=k, m=m, R=R[(k, m)], T=T, E_noise=En)
    return best, table


def run_l3(geo, vg, mode, sol0, G):
    rows = []
    for d in L3_DELTAS:
        I1, it1 = drain_current_at(geo, vg, mode, d, sol0, 1e-13)
        I2, it2 = drain_current_at(geo, vg, mode, d, sol0, 1e-15)
        ig4 = i_gca(geo, vg, d, 4); ig3 = i_gca(geo, vg, d, 3)
        r = I1 / ig4 - 1
        e_i = abs(I2 - I1) / abs(ig4) + abs(ig4 - ig3) / abs(ig4)
        rows.append(dict(delta_V=d, I2D_A_per_cm=I1, I2D_tight=I2, IGCA4=ig4, IGCA3=ig3, r=r, e_i=e_i,
                         newton_its=[it1, it2]))
    # back to V_D = 0
    ds.circuit_alter(name='V1', value=0.)
    if mode == 'a':
        set_lateral_data(geo, vg, sol0)
    D.converge()
    best, table = neville([x['delta_V'] for x in rows], [x['r'] for x in rows], [x['e_i'] for x in rows])
    out = dict(points=rows, tableau=table)
    if best is None:
        out.update(status='non-asymptotic', r0=None, E=None)
    else:
        out.update(status='ok', r0=best['R'], E=best['T'] + best['E_noise'], choice=best)
    return out


# ----------------------------------------------------------------------------- one (mode, VG)
def t1_point(geo, vg, mode, do_l3):
    r = dict(vg=vg, mode=mode)
    ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
    sol0, r['equilibrium'] = equilibrate(geo, vg, mode, TAU_EQ)
    G, ns = g_gca(geo, sol0)
    r['G_GCA_h_S_per_cm'] = G; r['Ns_FV_cm2'] = ns
    r['G_2D_midcolumn_S_per_cm'] = M.gca_conductance(geo.L_cm)['G_S_per_cm']
    eq = eqnums()
    gl1, gsec, nunk = l1(geo, eq); r['n_unknowns'] = nunk
    gns = newton_step_G()
    Y1, _ = l2(geo, F1)
    Y2, sec2 = l2(geo, F2)
    r.update(G_L1=gl1, G_L1_sections=gsec, G_NS=gns, Y_f1=[Y1.real, Y1.imag], Y_f2=[Y2.real, Y2.imag],
             ReY_sections_f2=[float(x.real) for x in sec2])
    r0_L1 = gl1 / G - 1; r0_L2 = Y2.real / G - 1
    # E_eq: tighten 2-D Newton and 1COL by 100
    rt = D.newton(TAU_EQ / 100, extra=0)
    sol_t = onecol_sol(geo, vg, tol=TAU_EQ / 100)
    Gt, _ = g_gca(geo, sol_t)
    gl1_t, _, _ = l1(geo, eqnums())
    Yt, _ = l2(geo, F2)
    E_eq_L1 = abs(gl1_t / Gt - 1 - r0_L1); E_eq_L2 = abs(Yt.real / Gt - 1 - r0_L2)
    E_freq = abs(Y1.real - Y2.real) / (0.99 * G)
    E_lin_L1 = max(abs(g - gl1) for g in gsec) / G + abs(gl1 - gns) / G
    E_lin_L2 = max(abs(float(x.real) - Y2.real) for x in sec2) / G + abs(Y2.real - gl1) / G
    E_L1 = E_eq_L1 + E_lin_L1; E_L2 = E_freq + E_eq_L2 + E_lin_L2
    r['tight'] = dict(newton_iterations=rt['iterations'], G_GCA_h=Gt, G_L1=gl1_t, ReY_f2=Yt.real)
    r['L1'] = dict(r0=r0_L1, E_eq=E_eq_L1, E_lin=E_lin_L1, E=E_L1)
    r['L2'] = dict(r0=r0_L2, E_freq=E_freq, E_eq=E_eq_L2, E_lin=E_lin_L2, E=E_L2)
    for k in ('L1', 'L2'):
        x = r[k]
        x['verdict'] = ('NON-RESOLVING (class I)' if x['E'] > 1e-6 else
                        'PASS' if abs(x['r0']) <= 3 * x['E'] else 'FAIL')
    if mode == 'a':
        ds.set_parameter(device=D.DEV, name='s_couple', value=1.)
        D.converge()
        gc, _, _ = l1(geo, eqnums())
        r['decoupling'] = dict(G_L1_s1=gc, rel_change=(gc - gl1) / G,
                               ok=abs(gc - gl1) / G <= E_L1)
        ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
        D.converge()
    if do_l3:
        r['L3'] = run_l3(geo, vg, mode, sol0, G)
        if r['L3']['status'] == 'ok':
            Ep = E_L2
            r['L3']['consistency_with_L2'] = dict(diff=abs(r0_L2 - r['L3']['r0']),
                                                  limit=3 * (Ep + r['L3']['E']),
                                                  ok=abs(r0_L2 - r['L3']['r0']) <= 3 * (Ep + r['L3']['E']))
    # boundary diagnostics
    A, *_ = export()
    r['continuity'] = continuity()
    r['gauss'] = gauss_closure(eq)
    r['oxide_lateral_single_entry_rows'] = oxide_lateral_rows(eq, A)['single_entry_rows']
    r['uniformity'] = [M.lateral_uniformity('Potential'), M.lateral_uniformity('Electrons')]
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--modes', default='a,b'); ap.add_argument('--vgs', default='0.8,1.25,2.0')
    ap.add_argument('--l3', type=int, default=0); ap.add_argument('--out', required=True)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    geo = Geo(a.L, a.level)
    rec = runlog.header(dict(L_um=a.L, level=a.level, modes=a.modes, vgs=a.vgs, l3=a.l3,
                             tau_eq=TAU_EQ, delta_L1=DELTA_L1, f1=F1, f2=F2, l3_deltas=L3_DELTAS,
                             pm=D.pm_checklist()))
    t0 = time.time()
    rec['mesh'] = D.build_mesh_explicit(geo.xs, geo.y_ox, geo.y_si, corner='tie')
    rec['mesh_quality'] = geo.quality()
    D.build_physics()
    P = replace(m0.M0Params.f0(), L=a.L * 1e-6)
    rec['points'] = {}
    for vg in [float(v) for v in a.vgs.split(',')]:
        for mode in a.modes.split(','):
            t = time.time()
            r = t1_point(geo, vg, mode, bool(a.l3))
            r['M0_G0_S_per_cm'] = m0.conductance0(vg, P) * 1e4
            r['wall_s'] = time.time() - t
            rec['points'][f'{mode}:VG={vg}'] = r
            print('T1', a.L, a.level, mode, vg, json.dumps(dict(L1=r['L1'], L2=r['L2'])), flush=True)
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(out, f't1_L{a.L:g}_lev{a.level}.json'), rec)
    return 0


if __name__ == '__main__':
    sys.exit(main())
