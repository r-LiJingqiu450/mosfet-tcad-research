"""C-6 exploratory probe W1: zero-bias conductance of a MOSFET channel bounded by n+ wells in
Ward's domain truncation ([WOC90] (2.4a), (2.5a,b); [W88] (5.2)), with an exact decomposition of
the source/drain end correction. EXPLORATORY (not the registered E-W1); F0 physics unchanged.

Boundary-value problem (DEVSIM coordinates: x along the channel, y = depth):
  * channel 0 < x < L, Si depth 0 < y < x_max = 1.2 um, oxide 20 nm, gate psi = V_G + psi_b0;
  * n+ well boundaries: Si edges x = 0 and x = L for 0 <= y <= h_w (well depth h_w, snapped to a
    depth line): psi = psi_bi = V_T ln(N_D/n_i) (N_D = 1e20) and phi_n = 0 / V_D (Dirichlet);
    the oxide triple points are tied to psi_bi (corner 'tie' construction of R-I(a));
  * lateral Si edges below h_w: homogeneous Neumann (psi and electrons);
  * bottom y = x_max: psi = psi_b0 over the full width; zero electron flux.
  The well interior is not meshed: Ward's equipotential n+ well with lambda_+/lambda >> 1.

Equilibrium: start from the exact laterally uniform 1COL state (lateral data = 1COL slice), then
homotopy of the lateral Dirichlet data t: (1-t) 1COL(y) + t psi_bi, t = 0 -> 1, adaptive steps,
Newton with the R05 D-3 rule (tau = 1e-12 V at t = 1).

Zero-bias conductance, three exact discrete quantities from the same equilibrium:
  At equilibrium the SG flux is J_e = q mu EIL V_T n_i B(v) e^{psi_0/V_T} (e^{-phi_1/V_T} - e^{-phi_0/V_T}),
  so dJ/dpsi = 0 and the linear response is a conductance network on the Si nodes with edge
  conductances c_e = q mu * EdgeCouple * EdgeInverseLength * n0 B(v), v = (psi0 - psi1)/V_T
  (symmetrized: (n0 B(v) + n1 B(-v))/2). Electrostatic feedback enters only at O(V_D^2).
  G_net   : network solve (Dirichlet phi_n on the well edges, natural elsewhere)    -- 2-D truth
  G_dev   : DEVSIM full coupled Jacobian linear response (equilibrated SuperLU, as T1 v2) -- cross-check
  G_slice : every column shorted (vertical conductances -> infinity): series sum of column
            conductances = Ward's quasi-1-D current assumption (6.33) with the EXACT 2-D charge.
            Rayleigh monotonicity: G_slice >= G_net (rigorous, discrete).
  G_rows  : vertical conductances -> 0: parallel 1-D rows.  G_rows <= G_net (rigorous).
  G_GCA   : q mu N_s / L from the 1COL slice on identical depth lines (= T1 reference).
  dL_X    = (L/2)(1 - G_GCA/G_X): end correction per end expressed as a length.
  dL_slice - dL_net >= 0 is the current-spreading part (Ward's (6.33) omits it by assumption).
Outputs also the surface potential psi(x, 0) and the column sheet density N(x).
"""
from __future__ import annotations

import argparse
import datetime as dt
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

from phase3 import meshmap as MM
from phase3.c1 import deck as D
from phase3.c3.t1 import onecol_sol, g_gca, _lut
from phase3.c4.t1v2 import systems, EqSolver, umf_solve, U_ROUND
from phase3.c1.corner_study import eqnums

ND = 1e20
PSI_BI = D.VT * math.log(ND / D.NI)
VERSION = 'wprobe-2026-10-31'


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds')


def lateral_map(L_um, a_um=0.35, h0_nm=0.5, hmax_nm=20.0):
    """Symmetric graded lateral map: geometric [0,a] fine at 0, uniform middle, mirrored [L-a,L]."""
    L, a, h0, hm = L_um * 1e-4, a_um * 1e-4, h0_nm * 1e-7, hmax_nm * 1e-7
    assert L > 2 * a
    best = None
    for n in range(10, 4000):
        Q = MM.stretch_for(a, h0, n)
        last = a * (Q - Q ** ((n - 1) / n)) / (Q - 1)
        if last <= hm:
            best = (n, Q); break
    n, Q = best
    nmid = max(1, math.ceil((L - 2 * a) / hm))
    return MM.AxisMap('lateral_w', (MM.Segment(0.0, a, n, 'geometric', Q),
                                    MM.Segment(a, L - a, nmid, 'uniform'),
                                    MM.Segment(L - a, L, n, 'geometric', 1.0 / Q)))


class GeoW:
    def __init__(self, L_um, level, h_w_um=0.5):
        self.L_um, self.level = L_um, level
        self.xmap = lateral_map(L_um)
        self.xs = self.xmap.nodes(level)
        self.L_cm = float(self.xs[-1])
        self.y_si = MM.f0_depth_map().nodes(level)
        self.y_ox = MM.oxide_map(20.0, 5).nodes(level)
        self.jw_si = int(np.argmin(np.abs(self.y_si - h_w_um * 1e-4)))   # snapped well depth
        self.h_w_cm = float(self.y_si[self.jw_si])

    def describe(self):
        return dict(L_um=self.L_um, level=self.level, n_lateral=int(len(self.xs)), n_depth_si=int(len(self.y_si)),
                    h_w_cm=self.h_w_cm, lateral_quality=self.xmap.quality(self.level),
                    dx_min_nm=float(np.diff(self.xs).min() * 1e7), dx_max_nm=float(np.diff(self.xs).max() * 1e7))


def build_mesh(geo: GeoW):
    xs = np.asarray(geo.xs, float); y_ox = np.asarray(geo.y_ox, float); y_si = np.asarray(geo.y_si, float)
    ys = np.concatenate([y_ox, y_si[1:]]); kif = len(y_ox) - 1
    nx, ny = len(xs), len(ys)
    nid = lambda i, j: j * nx + i
    coords = []
    for j in range(ny):
        for i in range(nx):
            coords += [float(xs[i]), float(ys[j]), 0.0]
    names = [D.OX, D.SI, 'ox_si', 'gate', 'source', 'drain', 'body', 'oxc_s', 'oxc_d']
    P = {n: k for k, n in enumerate(names)}
    el = []
    for j in range(ny - 1):
        reg = P[D.OX] if j < kif else P[D.SI]
        for i in range(nx - 1):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            el += [2, reg, a, b, c, 2, reg, a, c, d]
    for i in range(nx - 1):
        el += [1, P['ox_si'], nid(i, kif), nid(i + 1, kif)]
        el += [1, P['gate'], nid(i, 0), nid(i + 1, 0)]
    jw = kif + geo.jw_si
    for c, i in (('source', 0), ('drain', nx - 1)):
        for j in range(kif, jw):
            el += [1, P[c], nid(i, j), nid(i, j + 1)]
    for i in range(nx - 1):
        el += [1, P['body'], nid(i, ny - 1), nid(i + 1, ny - 1)]
    el += [1, P['oxc_s'], nid(0, kif), nid(0, kif)]
    el += [1, P['oxc_d'], nid(nx - 1, kif), nid(nx - 1, kif)]
    m = 'mw'
    ds.create_gmsh_mesh(mesh=m, coordinates=coords, elements=el, physical_names=names)
    ds.add_gmsh_region(gmsh_name=D.OX, mesh=m, region=D.OX, material='Ox')
    ds.add_gmsh_region(gmsh_name=D.SI, mesh=m, region=D.SI, material='Si')
    ds.add_gmsh_interface(gmsh_name='ox_si', mesh=m, name='ox_si', region0=D.OX, region1=D.SI)
    for c, r in (('gate', D.OX), ('source', D.SI), ('drain', D.SI), ('body', D.SI), ('oxc_s', D.OX), ('oxc_d', D.OX)):
        ds.add_gmsh_contact(gmsh_name=c, mesh=m, name=c, region=r, material='metal')
    ds.finalize_mesh(mesh=m)
    ds.create_device(mesh=m, device=D.DEV)
    D.CORNER['on'] = True; D.CORNER['mode'] = 'tie'
    return dict(nx=nx, ny=ny, kif=kif, jw=jw)


def contact_mask(region, geo, xs, ys):
    """Nodes whose PsiLat is the well value: Si well edges (y <= h_w) and the oxide triple points."""
    L = geo.L_cm
    edge = (np.abs(xs) < 1e-15) | (np.abs(xs - L) < 1e-15 * max(1, L))
    if region == D.SI:
        return edge & (ys <= geo.h_w_cm * (1 + 1e-12))
    return edge & (np.abs(ys) < 1e-15)


def set_psilat(t, geo, base):
    for r in (D.SI, D.OX):
        xs, ys = map(np.array, D.coords(r))
        m = contact_mask(r, geo, xs, ys)
        v = base[r].copy()
        v[m] = (1 - t) * base[r][m] + t * PSI_BI
        ds.set_node_values(device=D.DEV, region=r, name='PsiLat', values=list(v))


def equilibrium(geo, vg, tau=1e-12):
    sol0 = onecol_sol(geo, vg)
    lut = _lut(sol0)
    base = {}
    for r in (D.SI, D.OX):
        _, ys = D.coords(r)
        base[r] = np.array([lut[round(float(y), 15)] for y in ys])
        ds.set_node_values(device=D.DEV, region=r, name='PsiLat', values=list(base[r]))
        ds.set_node_values(device=D.DEV, region=r, name='Potential', values=list(base[r]))
    ds.set_node_values(device=D.DEV, region=D.SI, name='Electrons', values=[D.NI * math.exp(p / D.VT) for p in base[D.SI]])
    ds.set_parameter(device=D.DEV, name='gate_bias', value=vg)
    ds.circuit_alter(name='V1', value=0.)
    D.set_lateral_dirichlet(True)
    r0 = D.newton(1e-10)
    assert r0['converged'], r0['history'][-3:]
    log = [dict(t=0.0, it=r0['iterations'])]
    t, h = 0.0, 0.05
    while t < 1.0:
        tn = min(1.0, t + h)
        X = D._state()
        set_psilat(tn, geo, base)
        r = D.newton(1e-10, max_it=60)
        log.append(dict(t=tn, step=h, converged=r['converged'], it=r['iterations']))
        if r['converged']:
            t = tn; h = min(0.25, h * 1.5)
        else:
            D._restore(X); set_psilat(t, geo, base); h /= 2
            if h < 1e-4:
                raise D.NotConverged(f'homotopy failed at t={tn}: {log[-4:]}')
    rf = D.newton(tau, extra=2, max_it=60)
    assert rf['converged'], rf['history'][-3:]
    return sol0, dict(homotopy=log, final_iterations=rf['iterations'], final_updates=rf['history'][-3:])


def bern(v):
    v = np.asarray(v, float)
    out = np.empty_like(v)
    small = np.abs(v) < 1e-6
    out[small] = 1 - v[small] / 2 + v[small] ** 2 / 12
    vs = v[~small]
    out[~small] = vs / np.expm1(vs)
    return out


def network(geo):
    """Edge conductances of the zero-bias linear response (S per cm of width)."""
    ds.edge_from_node_model(device=D.DEV, region=D.SI, node_model='node_index')
    n0i = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='node_index@n0')).astype(int)
    n1i = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='node_index@n1')).astype(int)
    ec = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='EdgeCouple'))
    eil = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='EdgeInverseLength'))
    psi = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Potential'))
    n = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Electrons'))
    v = (psi[n0i] - psi[n1i]) / D.VT
    ne_a, ne_b = n[n0i] * bern(v), n[n1i] * bern(-v)
    c = D.Q * D.MU_N * ec * eil * 0.5 * (ne_a + ne_b)
    asym = float(np.max(np.abs(ne_a - ne_b) / np.maximum(np.abs(ne_a + ne_b) / 2, 1e-300)))
    return n0i, n1i, c, asym, psi, n


def solve_net(N, n0i, n1i, c, src, drn):
    """Dirichlet phi = 0 (src), 1 (drn); returns drain current (= G) and the source current."""
    W = sp.coo_matrix((np.concatenate([c, c]), (np.concatenate([n0i, n1i]), np.concatenate([n1i, n0i]))), shape=(N, N)).tocsr()
    deg = np.asarray(W.sum(axis=1)).ravel()
    Lap = (sp.diags(deg) - W).tocsr()
    fixed = np.zeros(N, bool); fixed[src] = True; fixed[drn] = True
    phi = np.zeros(N); phi[drn] = 1.0
    fr = np.where(~fixed)[0]
    A = Lap[fr][:, fr].tocsc(); b = -Lap[fr][:, np.where(fixed)[0]] @ phi[fixed]
    s = 1.0 / np.sqrt(A.diagonal())                       # symmetric Jacobi scaling
    As = (sp.diags(s) @ A @ sp.diags(s)).tocsc()
    y = spla.splu(As).solve(s * b)
    phi[fr] = s * y
    res = float(np.max(np.abs(A @ phi[fr] - b)) / np.max(np.abs(b)))
    I = Lap @ phi
    return float(np.sum(I[drn])), float(np.sum(I[src])), phi, res


def decomposition(geo, sol0):
    xs, ys = map(np.array, D.coords(D.SI))
    n0i, n1i, c, asym, psi, n = network(geo)
    N = len(xs)
    L = geo.L_cm
    src = np.where(contact_mask(D.SI, geo, xs, ys) & (xs < L / 2))[0]
    drn = np.where(contact_mask(D.SI, geo, xs, ys) & (xs > L / 2))[0]
    G_net, I_src, phi, res = solve_net(N, n0i, n1i, c, src, drn)
    # classify edges
    lat = np.abs(ys[n0i] - ys[n1i]) < 1e-18
    vert = np.abs(xs[n0i] - xs[n1i]) < 1e-18
    diag = ~(lat | vert)
    # G_slice: short every column -> series of column-to-column conductances
    xcol = np.array(sorted(set(np.round(xs, 16))))
    ci = np.searchsorted(xcol, np.round(xs, 16))
    a_, b_ = np.minimum(ci[n0i], ci[n1i]), np.maximum(ci[n0i], ci[n1i])
    between = (b_ - a_) == 1
    Ccol = np.bincount(a_[between], weights=c[between], minlength=len(xcol) - 1)[:len(xcol) - 1]
    G_slice = 1.0 / float(np.sum(1.0 / Ccol))
    # G_rows: cut vertical (and diagonal) edges; each row is a chain; only rows with both ends Dirichlet conduct
    yrow = np.array(sorted(set(np.round(ys, 18))))
    G_rows = 0.0
    for y in yrow:
        m = lat & (np.abs(ys[n0i] - y) < 1e-18)
        has_contact = np.any(np.abs(ys[src] - y) < 1e-18)
        if has_contact and m.any():
            G_rows += 1.0 / float(np.sum(1.0 / c[m]))
    G_GCA, Ns = g_gca(geo, sol0)
    # profiles
    jif = np.abs(ys) < 1e-18
    order = np.argsort(xs[jif])
    surf = dict(x_cm=xs[jif][order].tolist(), psi_V=psi[jif][order].tolist())
    w = dict(zip(np.round(sol0['y'][-len(geo.y_si):], 18), sol0['w_si']))
    wn = np.array([w.get(round(float(y), 18), np.nan) for y in ys])
    Ncol = np.bincount(ci, weights=n * wn, minlength=len(xcol))
    dl = lambda G: 0.5 * L * (1 - G_GCA / G) * 1e7       # nm per end
    return dict(G_GCA=G_GCA, Ns_GCA_cm2=Ns, G_net=G_net, G_net_source_current=I_src, net_rel_residual=res,
                G_slice=G_slice, G_rows=G_rows,
                dL_net_nm=dl(G_net), dL_slice_nm=dl(G_slice), dL_rows_nm=dl(G_rows),
                dL_spread_nm=dl(G_slice) - dl(G_net),
                bounds_ok=bool(G_rows <= G_net * (1 + 1e-12) and G_net <= G_slice * (1 + 1e-12)),
                kcl_rel=abs(G_net + I_src) / G_net, sg_equilibrium_asym=asym,
                n_edges=dict(lateral=int(lat.sum()), vertical=int(vert.sum()), diagonal=int(diag.sum()),
                             diagonal_max_coupling_rel=float(np.max(np.abs(c[diag])) / np.max(c)) if diag.any() else 0.0),
                surface=surf, column_sheet_density=dict(x_cm=xcol.tolist(), N_cm2=Ncol.tolist(), N_GCA=Ns))


def devsim_G():
    eq = eqnums()
    A, Ad, b0, bd = systems()
    n = A.shape[0]; iV = int(eq[('circuit', 'V1.I')][0])
    rhs = (bd - b0) / 1e-3
    s = EqSolver(A)
    x = s.solve(rhs)
    c = np.zeros(n); c[iV] = -1.0
    lam = s.solve_T(c)
    xm = umf_solve(A, rhs)
    r = A @ x - rhs
    G = -float(x[iV])
    sc = abs(A) @ np.abs(x) + np.abs(rhs)
    return dict(G_E=G, G_M=-float(xm[iV]), eta=s.eta(A, x, rhs), E_dwr=abs(float(lam @ r)) / abs(G),
                E_round=U_ROUND * float(np.linalg.norm(np.abs(lam) * sc)) / abs(G),
                kappa_G=float(np.abs(lam) @ sc) / abs(G), n_unknowns=int(n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--vg', type=float, required=True); ap.add_argument('--hw', type=float, default=0.5)
    ap.add_argument('--no-devsim', action='store_true',
                    help='skip the DEVSIM coupled-Jacobian cross-check (memory; level 2). The network G is the '
                         'primary quantity and was cross-validated against DEVSIM on all level-0/1 states')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f'w1_L{a.L:g}_lev{a.level}_VG{a.vg:g}_hw{a.hw:g}.json')
    rec = dict(version=VERSION, start_utc=now(), args=vars(a), psi_bi_V=PSI_BI, N_D=ND)
    t0 = time.time()
    geo = GeoW(a.L, a.level, a.hw)
    rec['geometry'] = geo.describe()
    rec['mesh'] = build_mesh(geo)
    D.build_physics()
    sol0, rec['equilibrium'] = equilibrium(geo, a.vg)
    rec['t_equilibrium_s'] = time.time() - t0
    rec['decomposition'] = decomposition(geo, sol0)
    dc = rec['decomposition']
    if a.no_devsim:
        rec['devsim'] = None
        rec['cross_check'] = dict(net_vs_devsim_rel=None, devsim_E_vs_M_rel=None, note='skipped (--no-devsim)')
        dv = {}
    else:
        rec['devsim'] = devsim_G()
        dv = rec['devsim']
        rec['cross_check'] = dict(net_vs_devsim_rel=(dc['G_net'] - dv['G_E']) / dv['G_E'],
                                  devsim_E_vs_M_rel=(dv['G_E'] - dv['G_M']) / dv['G_E'])
    rec['wall_s'] = time.time() - t0
    rec['end_utc'] = now()
    json.dump(rec, open(path, 'w'), indent=1)
    print(json.dumps({k: (v if not isinstance(v, (dict, list)) else '...') for k, v in dc.items() if k not in ('surface', 'column_sheet_density')}, indent=1))
    print('devsim', json.dumps(dv), 'cross', json.dumps(rec['cross_check']), 'wall', rec['wall_s'])


if __name__ == '__main__':
    main()
