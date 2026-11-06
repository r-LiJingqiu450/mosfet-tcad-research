"""Shared deck-R construction for WP2 paths P1 (network) and P2 (drift-diffusion).

Device definition: benchmark/specs/WP3_Preregistration_AR_and_HeldOut_2026-11-04.md, section 4.
  * resolved n+ wells, ohmic far-edge contacts (psi = V_T asinh(N_net/2n_i) + s_couple*V, phi_n = V), s_couple = 1;
  * donors N_D(x, y) = 1e20 f_x(x) f_y(y): lateral erfc (metallurgical junction at x = 0, L) and vertical erfc
    centred at h_w = 0.5 um (sigma = sigma_y = 15 nm)  [pre-registration star item];
  * body at y = 1.2 um: psi = psi_b0 and either B0 (zero electron flux) or B1 (ohmic: n = n_i exp(psi/V_T));
  * meshes: phase3/c7/endprobe.Geo deck-R maps (rprobe.lateral_nodes, endprobe.depth_map, oxide_map(t_ox, 5)).
    The contact extent is the depth line nearest to h_w (as in the existing deck R); its effect is bounded
    separately (benchmark/wp2/contact_extent_check.py).
Reuses the existing deck (phase3/c1/deck.py) unchanged; nothing in phase3/ is modified.
"""
from __future__ import annotations

import math
import os
import sys
from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import erfc, erfcinv

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import devsim as ds  # noqa: E402

from phase3.c1 import deck as D  # noqa: E402
from phase3.c3.t1 import onecol_sol, g_gca, _lut  # noqa: E402
from phase3.c6 import wprobe as WP  # noqa: E402
from phase3.c7 import endprobe as EP  # noqa: E402

ND = 1e20
HW_UM, LW_UM, SIGMA_NM, SIGMA_Y_NM, XMAX_UM = 0.5, 0.3, 15.0, 15.0, 1.2

DEVICES = {
    'dev1': dict(na=1e18, tox=2.0),
    'dev2': dict(na=1e16, tox=20.0),
}


def make_geo(L, level, tox, na, vg):
    a = SimpleNamespace(deck='R', L=L, level=level, tox=tox, na=na, ov=0.0, sigma=SIGMA_NM, lw=LW_UM, hw=HW_UM,
                        vg=vg, psibi_shift=0.0)
    EP.set_bulk(na)
    geo = EP.Geo(a)
    assert abs(geo.y_si[-1] - XMAX_UM * 1e-4) < 1e-15
    return a, geo


def donor(geo, xs, ys):
    """N_D(x, y) = 1e20 f_x(x) f_y(y) (cm^-3), pre-registration section 4."""
    s2 = math.sqrt(2.0) * SIGMA_NM * 1e-7
    xc = -s2 * float(erfcinv(2.0 * geo.na / ND))
    L = geo.L_cm
    fx = 0.5 * erfc((xs - xc) / s2) + 0.5 * erfc(((L - xs) - xc) / s2)
    fy = 0.5 * erfc((ys - HW_UM * 1e-4) / (math.sqrt(2.0) * SIGMA_Y_NM * 1e-7))
    return ND * fx * fy


def terminal_masks(geo, xs, ys):
    x0, x1, yb = float(np.min(xs)), float(np.max(xs)), float(np.max(ys))
    tol = 1e-15
    in_contact = ys <= geo.h_w_cm * (1 + 1e-12)
    return dict(source=(np.abs(xs - x0) < tol) & in_contact,
                drain=(np.abs(xs - x1) < tol) & in_contact,
                body=np.abs(ys - yb) < tol)


def add_body_electron_contact():
    """B1: ohmic body for electrons, phi_n = 0  ->  n = n_i exp(Potential/V_t) on the body contact."""
    eq = 'Electrons - n_i*exp(Potential/V_t)'
    ds.contact_node_model(device=D.DEV, contact='body', name='body_n', equation=eq)
    ds.contact_node_model(device=D.DEV, contact='body', name='body_n:Electrons', equation='1')
    ds.contact_node_model(device=D.DEV, contact='body', name='body_n:Potential', equation='-n_i*exp(Potential/V_t)/V_t')
    ds.contact_equation(device=D.DEV, contact='body', name='ElectronContinuityEquation', node_model='body_n',
                        edge_current_model='ElectronCurrent')


def build(L, level, tox, na, vg, body):
    """Mesh + physics + equilibrium. Returns (args, geo, sol0, record)."""
    a, geo = make_geo(L, level, tox, na, vg)
    rec = dict(bulk=dict(N_A=D.NA, P0=D.P0, psi_b0=D.PSI_B0))
    rec['geometry'] = geo.describe()
    rec['geometry']['contact_extent_cm'] = geo.h_w_cm
    rec['geometry']['h_w_doping_cm'] = HW_UM * 1e-4
    rec['mesh'] = EP.build_mesh(geo)
    D.build_physics()
    ds.set_parameter(device=D.DEV, name='s_couple', value=1.0)
    if body == 'B1':
        add_body_electron_contact()
    sol0, eq = equilibrium(geo, a)
    rec['equilibrium'] = eq
    k = sol0['k']
    rec['onecol'] = dict(psi_s_V=float(sol0['psi'][k]), psi_hat=float((sol0['psi'][k] - D.PSI_B0) / D.VT),
                         converged=bool(sol0['converged']), last_update_V=float(sol0['last_update_V']))
    return a, geo, sol0, rec


def equilibrium(geo, a, tau=1e-12):
    """Deck-R equilibrium as phase3/c7/endprobe.equilibrium (deck R branch), with the WP2 donor profile."""
    sol0 = onecol_sol(geo, a.vg)
    lut = _lut(sol0)
    xs, ys = map(np.array, D.coords(D.SI)); xo, yo = map(np.array, D.coords(D.OX))
    psi_si = np.array([lut[round(float(y), 15)] for y in ys])
    psi_ox = np.array([lut[round(float(y), 15)] for y in yo])
    nd = donor(geo, xs, ys)
    net = nd - geo.na
    ds.node_solution(device=D.DEV, region=D.SI, name='NetDop')
    ds.set_node_values(device=D.DEV, region=D.SI, name='NetDop', values=list(net))
    D._model(D.SI, 'PotentialNodeCharge', '-ElectronCharge*kahan3(Holes, -Electrons, NetDop)', ('Potential', 'Electrons'))
    neutral = D.VT * np.arcsinh(net / (2 * D.NI))
    ntype = net > 0
    psi_si[ntype] = neutral[ntype]
    T = terminal_masks(geo, xs, ys)
    cs = T['source'] | T['drain']
    psi_w = float(np.max(neutral[cs]))
    over_well = (xo < -1e-15) | (xo > geo.L_cm * (1 + 1e-12))
    psi_ox[over_well & (np.abs(yo) < 1e-15)] = psi_w
    ox_corner = ((np.abs(xo - xs.min()) < 1e-15) | (np.abs(xo - xs.max()) < 1e-15)) & (np.abs(yo) < 1e-15)
    for r, v, m, lat in ((D.SI, psi_si, cs, neutral), (D.OX, psi_ox, ox_corner, np.full(len(psi_ox), psi_w))):
        ds.set_node_values(device=D.DEV, region=r, name='Potential', values=list(v))
        ds.set_node_values(device=D.DEV, region=r, name='PsiLat', values=list(np.where(m, lat, v)))
    ds.set_node_values(device=D.DEV, region=D.SI, name='Electrons', values=list(D.NI * np.exp(psi_si / D.VT)))
    ds.set_parameter(device=D.DEV, name='gate_bias', value=a.vg)
    ds.circuit_alter(name='V1', value=0.)
    D.set_lateral_dirichlet(True)
    r = D.newton(tau, max_it=100, extra=2)
    log = dict(n_donor_nodes=int(np.sum(nd > geo.na)), psi_contact_V=psi_w, tau_V=tau,
               converged=bool(r['converged']), iterations=int(r['iterations']), last=r['history'][-3:])
    if not r['converged']:
        raise D.NotConverged(f'equilibrium failed: {r["history"][-4:]}')
    return sol0, log


# --------------------------------------------------------------------------------------------- P3
def eqnums_capture():
    """Equation numbering of the current Newton system without factoring the raw matrix.

    phase3/c1/corner_study.eqnums captures the numbering inside a solver callback that also factors the RAW
    (unscaled) Jacobian with SuperLU; on dev2 level 1 that factorisation reports 'exactly singular' although
    DEVSIM's own solver and the equilibrated solve handle the system. Only the numbering is needed, so this
    callback returns a zero update; the state is restored afterwards exactly as in corner_study.eqnums."""
    from phase3.c1.c1b_jacobian import equation_numbers, snapshot, restore
    rec = {}

    def cb(**kw):
        act = kw['action']
        if act == 'init':
            rec['n'] = kw['n']
            return dict(solver_object=object(), matrix_format='csc', status=True, message='')
        if act == 'factor':
            rec['eqnums'] = equation_numbers()
            return dict(status=True, message='')
        if act == 'solve':
            return dict(status=True, message='', x=[0.0] * rec['n'])
        raise RuntimeError(act)

    stock_cb = ds.get_parameter(name='solver_callback')
    try:
        stock_solver = ds.get_parameter(name='direct_solver')
    except ds.error:
        stock_solver = None
    X0 = snapshot()
    ds.set_parameter(name='direct_solver', value='custom')
    ds.set_parameter(name='solver_callback', value=cb)
    try:
        D.solve(abs_err=1e300, rel_err=1e300, it=1)
    finally:
        ds.set_parameter(name='solver_callback', value=stock_cb)
        if stock_solver is not None:
            ds.set_parameter(name='direct_solver', value=stock_solver)
    restore(X0)
    return rec['eqnums']


def devsim_smallsignal():
    """Same computation as phase3/c6/wprobe.devsim_G (G_E: equilibrated SuperLU; G_M: UMFPACK; DWR, round-off and
    kappa_G diagnostics), with the numbering from eqnums_capture()."""
    from phase3.c4.t1v2 import systems, EqSolver, umf_solve, U_ROUND
    eq = eqnums_capture()
    A, Ad, b0, bd = systems()
    n = A.shape[0]; iV = int(eq[('circuit', 'V1.I')][0])
    rhs = (bd - b0) / 1e-3
    s = EqSolver(A)
    x = s.solve(rhs)
    cvec = np.zeros(n); cvec[iV] = -1.0
    lam = s.solve_T(cvec)
    xm = umf_solve(A, rhs)
    r = A @ x - rhs
    G = -float(x[iV])
    sc = abs(A) @ np.abs(x) + np.abs(rhs)
    return dict(G_E=G, G_M=-float(xm[iV]), eta=s.eta(A, x, rhs), E_dwr=abs(float(lam @ r)) / abs(G),
                E_round=U_ROUND * float(np.linalg.norm(np.abs(lam) * sc)) / abs(G),
                kappa_G=float(np.abs(lam) @ sc) / abs(G), n_unknowns=int(n))


# --------------------------------------------------------------------------------------------- network
def laplacian(N, n0i, n1i, c):
    W = sp.coo_matrix((np.concatenate([c, c]), (np.concatenate([n0i, n1i]), np.concatenate([n1i, n0i]))),
                      shape=(N, N)).tocsr()
    return (sp.diags(np.asarray(W.sum(axis=1)).ravel()) - W).tocsr()


def terminal_solve(Lap, n0i, n1i, c, terminals, excite):
    """phi = 1 on terminal `excite`, 0 on the other terminals; natural (zero-flux) elsewhere.
    Returns terminal currents I_t = sum_{i in t} (L phi)_i, the edge energy phi^T L phi and the solution."""
    N = Lap.shape[0]
    fixed = np.zeros(N, bool)
    for idx in terminals.values():
        fixed[idx] = True
    phi = np.zeros(N); phi[terminals[excite]] = 1.0
    fr = np.where(~fixed)[0]
    A = Lap[fr][:, fr].tocsc(); b = -Lap[fr][:, np.where(fixed)[0]] @ phi[fixed]
    s = 1.0 / np.sqrt(A.diagonal())
    As = (sp.diags(s) @ A @ sp.diags(s)).tocsc()
    phi[fr] = s * spla.splu(As).solve(s * b)
    res = float(np.max(np.abs(A @ phi[fr] - b)) / np.max(np.abs(b)))
    I = Lap @ phi
    out = {t: float(np.sum(I[idx])) for t, idx in terminals.items()}
    energy = float(np.sum(c * (phi[n0i] - phi[n1i]) ** 2))
    return dict(I=out, energy=energy, residual=res, phi=phi)
