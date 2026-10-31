"""C-7 end-correction probe: zero-bias source/drain end correction dL(V_G) of a MOSFET channel,
generalising phase3/c6/{wprobe,rprobe}.py (regression: identical results for the C-6 settings).

Decks (DEVSIM x along the channel, y = depth; F0 physics: Boltzmann, frozen holes, constant mu_n,
no G-R; N_A, t_ox selectable):
  W  Ward truncation: channel 0<x<L; well edges x=0, L for 0<=y<=h_w carry psi = psi_bi + psi_shift,
     phi_n Dirichlet; oxide triple points tied; Neumann below h_w; bottom psi = psi_b0.
     psi_shift emulates the built-in-potential changes of degenerate statistics / band-gap narrowing
     in the n+ well (sensitivity only).
  R  resolved n+ wells: Si x in [-l_w, L+l_w]; N_D(x, y) = N_D * lateral profile * [y <= h_w];
     abrupt (sigma = 0, box) or lateral erfc grading with standard deviation sigma, positioned so
     that the METALLURGICAL junction (N_D = N_A) is at x = 0 and x = L; gate on the oxide top over
     [-ov, L + ov] (overlap ov >= 0); ohmic far-edge contacts psi = V_T asinh(N_net/2n_i), phi_n.
Zero-bias conductance (C-6 P4: exact network with equilibrium edge conductances):
  constant mobility: G_net, G_slice (columns shorted, upper bound), G_rows (lower bound), 1COL G_GCA.
  field-dependent mobility (post-processing, exact under P4 because d(mu)/dV multiplies grad(phi_n)=0):
     Lombardi-type surface mobility 1/mu = 1/mu_b + 1/mu_ac + 1/mu_sr with
     mu_ac = B/E_perp + C N^lam /(T E_perp^(1/3)), mu_sr = delta/E_perp^2 (B=4.75e7, C=1.74e5, lam=0.125,
     delta=5.82e14, T=300; no distance damping); E_perp = |E_y| for lateral edges, |E_x| for vertical
     edges (node averages of edge fields); mu_b = 1350 (F0).  GCA reference for the mobility case:
     same-mesh sheet conductance of the channel-centre column, G_GCA,mesh = sigma_mid / L.
  dL_X = (L/2)(1 - G_GCA/G_X) per end, L = metallurgical channel length.
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
from scipy.special import erfc, erfcinv

import devsim as ds

from phase3 import meshmap as MM
from phase3.c1 import deck as D
from phase3.c3.t1 import onecol_sol, g_gca, _lut
from phase3.c6 import wprobe as WP

VERSION = 'endprobe-2026-10-31'
ND = 1e20
LOMB = dict(B=4.75e7, C=1.74e5, lam=0.125, delta=5.82e14, T=300.0, mu_b=1350.0)


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds')


def set_bulk(na):
    D.NA = na
    D.P0 = .5 * (na + math.sqrt(na * na + 4 * D.NI * D.NI))
    D.PSI_B0 = -D.VT * math.log(D.P0 / D.NI)


def depth_map(na, h0_nm=0.25, hmax_nm=20.0, xmax_um=1.2):
    """F0 depth map with breakpoint 3 W_dep(N_A) (identical to MM.f0_depth_map for N_A = 1e16)."""
    xm, xb = xmax_um * 1e-4, 3 * MM.wdep_max_cm(na_cm3=na)
    h0, hm = h0_nm * 1e-7, hmax_nm * 1e-7
    for n in range(20, 2000):
        Q = MM.stretch_for(xb, h0, n)
        if xb * (Q - Q ** ((n - 1) / n)) / (Q - 1) <= hm:
            break
    return MM.AxisMap('si_depth', (MM.Segment(0.0, xb, n, 'geometric', Q),
                                   MM.Segment(xb, xm, math.ceil((xm - xb) / hm), 'uniform')))


class Geo:
    def __init__(self, a):
        self.deck, self.L_um, self.level = a.deck, a.L, a.level
        self.tox_nm, self.na, self.ov_cm, self.sigma_cm = a.tox, a.na, a.ov * 1e-7, a.sigma * 1e-7
        self.lw_cm = a.lw * 1e-4
        core = WP.lateral_map(a.L).nodes(a.level)
        L = float(core[-1])
        if a.deck == 'W':
            xs = core; self.i0 = 0
        else:
            from phase3.c6.rprobe import lateral_nodes
            xs, self.i0 = lateral_nodes(a.L, a.level, a.lw)
            if self.ov_cm > 0:
                xs = np.unique(np.concatenate([xs, [-self.ov_cm, L + self.ov_cm]]))
                self.i0 = int(np.argmin(np.abs(xs)))
        self.xs = xs
        self.i1 = int(np.argmin(np.abs(xs - L)))
        assert xs[self.i0] == 0.0 and abs(xs[self.i1] - L) < 1e-18
        self.L_cm = L
        self.y_si = depth_map(a.na).nodes(a.level)
        self.y_ox = MM.oxide_map(a.tox, 5).nodes(a.level)
        self.jw_si = int(np.argmin(np.abs(self.y_si - a.hw * 1e-4)))
        self.h_w_cm = float(self.y_si[self.jw_si])

    def describe(self):
        h = np.diff(self.xs)
        return dict(deck=self.deck, L_um=self.L_um, level=self.level, t_ox_nm=self.tox_nm, N_A=self.na,
                    ov_cm=self.ov_cm, sigma_cm=self.sigma_cm, l_w_cm=self.lw_cm, h_w_cm=self.h_w_cm,
                    n_lateral=int(len(self.xs)), n_depth_si=int(len(self.y_si)), n_ox=int(len(self.y_ox)),
                    dx_min_nm=float(h.min() * 1e7), dx_max_nm=float(h.max() * 1e7),
                    max_adjacent_ratio=float(np.max(np.maximum(h[1:] / h[:-1], h[:-1] / h[1:]))))


def build_mesh(geo):
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
    if geo.deck == 'W':
        g0, g1 = 0, nx - 1
    else:
        g0 = int(np.argmin(np.abs(xs + geo.ov_cm))); g1 = int(np.argmin(np.abs(xs - geo.L_cm - geo.ov_cm)))
    for i in range(nx - 1):
        el += [1, P['ox_si'], nid(i, kif), nid(i + 1, kif)]
    for i in range(g0, g1):
        el += [1, P['gate'], nid(i, 0), nid(i + 1, 0)]
    jw = kif + geo.jw_si
    for c, i in (('source', 0), ('drain', nx - 1)):
        for j in range(kif, jw):
            el += [1, P[c], nid(i, j), nid(i, j + 1)]
    for i in range(nx - 1):
        el += [1, P['body'], nid(i, ny - 1), nid(i + 1, ny - 1)]
    el += [1, P['oxc_s'], nid(0, kif), nid(0, kif)]
    el += [1, P['oxc_d'], nid(nx - 1, kif), nid(nx - 1, kif)]
    m = 'm7'
    ds.create_gmsh_mesh(mesh=m, coordinates=coords, elements=el, physical_names=names)
    ds.add_gmsh_region(gmsh_name=D.OX, mesh=m, region=D.OX, material='Ox')
    ds.add_gmsh_region(gmsh_name=D.SI, mesh=m, region=D.SI, material='Si')
    ds.add_gmsh_interface(gmsh_name='ox_si', mesh=m, name='ox_si', region0=D.OX, region1=D.SI)
    for c, r in (('gate', D.OX), ('source', D.SI), ('drain', D.SI), ('body', D.SI), ('oxc_s', D.OX), ('oxc_d', D.OX)):
        ds.add_gmsh_contact(gmsh_name=c, mesh=m, name=c, region=r, material='metal')
    ds.finalize_mesh(mesh=m)
    ds.create_device(mesh=m, device=D.DEV)
    D.CORNER['on'] = True; D.CORNER['mode'] = 'tie'
    return dict(nx=nx, ny=ny, kif=kif, jw=jw, gate_cols=[g0, g1])


def donor(geo, xs, ys):
    """Donor density (cm^-3) of the resolved wells."""
    inbox = ys <= geo.h_w_cm * (1 + 1e-12)
    L = geo.L_cm
    if geo.sigma_cm == 0.0:
        lat = ((xs < -1e-15) | (xs > L * (1 + 1e-12))).astype(float)
        return ND * lat * inbox
    s2 = math.sqrt(2.0) * geo.sigma_cm
    xc = -s2 * float(erfcinv(2.0 * geo.na / ND))          # profile centre, metallurgical junction at 0
    left = 0.5 * erfc((xs - xc) / s2); right = 0.5 * erfc(((L - xs) - xc) / s2)
    return ND * (left + right) * inbox


def contact_mask(region, geo, xs, ys):
    x0, x1 = float(geo.xs[0]), float(geo.xs[-1])
    edge = (np.abs(xs - x0) < 1e-15) | (np.abs(xs - x1) < 1e-15)
    if region == D.SI:
        return edge & (ys <= geo.h_w_cm * (1 + 1e-12))
    return edge & (np.abs(ys) < 1e-15)


def equilibrium(geo, a, tau=1e-12):
    sol0 = onecol_sol(geo, a.vg)
    lut = _lut(sol0)
    xs, ys = map(np.array, D.coords(D.SI)); xo, yo = map(np.array, D.coords(D.OX))
    psi_si = np.array([lut[round(float(y), 15)] for y in ys])
    psi_ox = np.array([lut[round(float(y), 15)] for y in yo])
    log = dict()
    if geo.deck == 'W':
        target = WP.PSI_BI + a.psibi_shift
        base = {D.SI: psi_si, D.OX: psi_ox}
        for r in (D.SI, D.OX):
            ds.set_node_values(device=D.DEV, region=r, name='PsiLat', values=list(base[r]))
            ds.set_node_values(device=D.DEV, region=r, name='Potential', values=list(base[r]))
        ds.set_node_values(device=D.DEV, region=D.SI, name='Electrons', values=list(D.NI * np.exp(psi_si / D.VT)))
        ds.set_parameter(device=D.DEV, name='gate_bias', value=a.vg)
        ds.circuit_alter(name='V1', value=0.)
        D.set_lateral_dirichlet(True)
        r0 = D.newton(1e-10)
        assert r0['converged'], r0['history'][-3:]
        hist = [dict(t=0.0, it=r0['iterations'])]
        t, h = 0.0, 0.05
        while t < 1.0:
            tn = min(1.0, t + h)
            X = D._state()
            for r in (D.SI, D.OX):
                xr, yr = map(np.array, D.coords(r)); m = contact_mask(r, geo, xr, yr)
                v = base[r].copy(); v[m] = (1 - tn) * base[r][m] + tn * target
                ds.set_node_values(device=D.DEV, region=r, name='PsiLat', values=list(v))
            r = D.newton(1e-10, max_it=60)
            hist.append(dict(t=tn, step=h, converged=r['converged'], it=r['iterations']))
            if r['converged']:
                t = tn; h = min(0.25, h * 1.5)
            else:
                D._restore(X); h /= 2
                for rr in (D.SI, D.OX):
                    xr, yr = map(np.array, D.coords(rr)); m = contact_mask(rr, geo, xr, yr)
                    v = base[rr].copy(); v[m] = (1 - t) * base[rr][m] + t * target
                    ds.set_node_values(device=D.DEV, region=rr, name='PsiLat', values=list(v))
                if h < 1e-4:
                    raise D.NotConverged(f'homotopy failed at t={tn}')
        log['homotopy'] = hist
    else:
        nd = donor(geo, xs, ys)
        net = nd - geo.na
        ds.node_solution(device=D.DEV, region=D.SI, name='NetDop')
        ds.set_node_values(device=D.DEV, region=D.SI, name='NetDop', values=list(net))
        D._model(D.SI, 'PotentialNodeCharge', '-ElectronCharge*kahan3(Holes, -Electrons, NetDop)', ('Potential', 'Electrons'))
        neutral = D.VT * np.arcsinh(net / (2 * D.NI))
        ntype = net > 0
        psi_si[ntype] = neutral[ntype]
        cs = contact_mask(D.SI, geo, xs, ys)
        psi_w = float(np.max(neutral[cs]))
        over_well = (xo < -1e-15) | (xo > geo.L_cm * (1 + 1e-12))
        psi_ox[over_well & (np.abs(yo) < 1e-15)] = psi_w
        for r, v, m in ((D.SI, psi_si, cs), (D.OX, psi_ox, contact_mask(D.OX, geo, xo, yo))):
            ds.set_node_values(device=D.DEV, region=r, name='Potential', values=list(v))
            pl = np.where(m, (neutral if r == D.SI else np.full(len(v), psi_w)), v)
            ds.set_node_values(device=D.DEV, region=r, name='PsiLat', values=list(pl))
        ds.set_node_values(device=D.DEV, region=D.SI, name='Electrons', values=list(D.NI * np.exp(psi_si / D.VT)))
        ds.set_parameter(device=D.DEV, name='gate_bias', value=a.vg)
        ds.circuit_alter(name='V1', value=0.)
        D.set_lateral_dirichlet(True)
        log['n_donor_nodes'] = int(np.sum(nd > geo.na))
        log['psi_contact_V'] = psi_w
    r = D.newton(tau, max_it=100, extra=2)
    log['final'] = dict(converged=r['converged'], iterations=r['iterations'], last=r['history'][-3:])
    if not r['converged']:
        raise D.NotConverged(f'equilibrium failed: {r["history"][-4:]}')
    return sol0, log


def node_fields(xs, ys, psi, n0i, n1i):
    """Node |E_y| and |E_x| from edge differences (tensor mesh; diagonals ignored)."""
    N = len(xs)
    lat = np.abs(ys[n0i] - ys[n1i]) < 1e-18
    ver = np.abs(xs[n0i] - xs[n1i]) < 1e-18
    out = {}
    for name, m, coord in (('Ex', lat, xs), ('Ey', ver, ys)):
        e = -(psi[n1i[m]] - psi[n0i[m]]) / (coord[n1i[m]] - coord[n0i[m]])
        s = np.bincount(n0i[m], weights=e, minlength=N) + np.bincount(n1i[m], weights=e, minlength=N)
        c = np.bincount(n0i[m], minlength=N) + np.bincount(n1i[m], minlength=N)
        out[name] = np.abs(s / np.maximum(c, 1))
    return out, lat, ver


def lombardi(E, Ndop):
    E = np.maximum(E, 1.0)
    p = LOMB
    mu_ac = p['B'] / E + p['C'] * np.power(np.maximum(Ndop, 1.0), p['lam']) / (p['T'] * np.cbrt(E))
    mu_sr = p['delta'] / E ** 2
    return 1.0 / (1.0 / p['mu_b'] + 1.0 / mu_ac + 1.0 / mu_sr)


def conductances(geo, c, n0i, n1i, xs, ys, L):
    """G_net, G_slice, G_rows and the same-mesh centre-column sheet conductance for edge conductances c."""
    src = np.where(contact_mask(D.SI, geo, xs, ys) & (xs < L / 2))[0]
    drn = np.where(contact_mask(D.SI, geo, xs, ys) & (xs > L / 2))[0]
    G, Is, phi, res = WP.solve_net(len(xs), n0i, n1i, c, src, drn)
    xcol = np.array(sorted(set(np.round(xs, 16))))
    ci = np.searchsorted(xcol, np.round(xs, 16))
    a_, b_ = np.minimum(ci[n0i], ci[n1i]), np.maximum(ci[n0i], ci[n1i])
    btw = (b_ - a_) == 1
    Ccol = np.bincount(a_[btw], weights=c[btw], minlength=len(xcol) - 1)[:len(xcol) - 1]
    G_slice = 1.0 / float(np.sum(1.0 / Ccol))
    lat = np.abs(ys[n0i] - ys[n1i]) < 1e-18
    G_rows = 0.0
    for y in np.array(sorted(set(np.round(ys[src], 18)))):
        m = lat & (np.abs(ys[n0i] - y) < 1e-18)
        if m.any():
            G_rows += 1.0 / float(np.sum(1.0 / c[m]))
    k = int(np.argmin(np.abs(0.5 * (xcol[:-1] + xcol[1:]) - L / 2)))
    sigma_mid = float(Ccol[k] * (xcol[k + 1] - xcol[k]))
    return dict(G_net=G, G_slice=G_slice, G_rows=G_rows, src_current=Is, net_rel_residual=res,
                sigma_mid=sigma_mid, bounds_ok=bool(G_rows <= G * (1 + 1e-12) and G <= G_slice * (1 + 1e-12)))


def analyse(geo, sol0):
    xs, ys = map(np.array, D.coords(D.SI))
    n0i, n1i, c, asym, psi, n = WP.network(geo)
    L = geo.L_cm
    G_GCA, Ns = g_gca(geo, sol0)
    dl = lambda Gref, G: 0.5 * L * (1 - Gref / G) * 1e7
    out = dict(G_GCA_1COL=G_GCA, Ns_cm2=Ns, sg_equilibrium_asym=asym)
    cst = conductances(geo, c, n0i, n1i, xs, ys, L)
    G_GCA_mesh = cst['sigma_mid'] / L
    cst.update(G_GCA_mesh_centre=G_GCA_mesh, centre_vs_1col=G_GCA_mesh / G_GCA - 1,
               dL_net_nm=dl(G_GCA, cst['G_net']), dL_slice_nm=dl(G_GCA, cst['G_slice']), dL_rows_nm=dl(G_GCA, cst['G_rows']),
               dL_net_meshref_nm=dl(G_GCA_mesh, cst['G_net']))
    cst['dL_spread_nm'] = cst['dL_slice_nm'] - cst['dL_net_nm']
    # diagnostic (C-7): slice conductance using only the surface layer y <= y_c (not a bound): tests whether the
    # slice-minus-network difference is carried by charge below the surface layer (junction spill-over)
    lat_e = np.abs(ys[n0i] - ys[n1i]) < 1e-18
    xcol = np.array(sorted(set(np.round(xs, 16)))); ci = np.searchsorted(xcol, np.round(xs, 16))
    a_, b_ = np.minimum(ci[n0i], ci[n1i]), np.maximum(ci[n0i], ci[n1i])
    for yc in (20e-7, 50e-7):
        m = lat_e & ((b_ - a_) == 1) & (ys[n0i] <= yc * (1 + 1e-12))
        Cc = np.bincount(a_[m], weights=c[m], minlength=len(xcol) - 1)[:len(xcol) - 1]
        k = int(np.argmin(np.abs(0.5 * (xcol[:-1] + xcol[1:]) - L / 2)))
        Gt = 1.0 / float(np.sum(1.0 / Cc)); Gref_t = Cc[k] * (xcol[k + 1] - xcol[k]) / L
        cst[f'dL_slice_top{int(round(yc * 1e7))}nm_nm'] = dl(Gref_t, Gt)
    out['const_mu'] = cst
    # field-dependent mobility
    F, lat, ver = node_fields(xs, ys, psi, n0i, n1i)
    nd = donor(geo, xs, ys) if geo.deck == 'R' else np.zeros(len(xs))
    Ndop = geo.na + nd
    Eperp = np.where(lat, 0.5 * (F['Ey'][n0i] + F['Ey'][n1i]), 0.5 * (F['Ex'][n0i] + F['Ex'][n1i]))
    mu = lombardi(Eperp, 0.5 * (Ndop[n0i] + Ndop[n1i]))
    cm = c * mu / D.MU_N
    mob = conductances(geo, cm, n0i, n1i, xs, ys, L)
    # far-field (1COL) mobility-weighted sheet conductance on the same depth lines and the same node-field rule
    ysi = np.asarray(sol0['y'][-len(geo.y_si):], float); psi1 = np.asarray(sol0['psi'][-len(geo.y_si):], float)
    n1 = np.asarray(sol0['n'][-len(geo.y_si):], float); w1 = np.asarray(sol0['w_si'][-len(geo.y_si):], float)
    ef = -np.diff(psi1) / np.diff(ysi)
    Ey1 = np.abs(np.concatenate([[ef[0]], 0.5 * (ef[:-1] + ef[1:]), [ef[-1]]]))
    mu1 = lombardi(Ey1, np.full(len(Ey1), geo.na))
    Gm = D.Q * float(np.sum(mu1 * n1 * w1)) / L
    mob['G_GCA_mesh_centre'] = mob['sigma_mid'] / L
    mob.update(G_GCA_1COL_mu=Gm, dL_net_nm=dl(Gm, mob['G_net']), dL_slice_nm=dl(Gm, mob['G_slice']),
               dL_rows_nm=dl(Gm, mob['G_rows']), mu_min=float(mu.min()), mu_max=float(mu.max()))
    mob['dL_spread_nm'] = mob['dL_slice_nm'] - mob['dL_net_nm']
    jif = np.abs(ys) < 1e-18
    o = np.argsort(xs[jif])
    out['lombardi_mu'] = mob
    out['surface'] = dict(x_cm=xs[jif][o].tolist(), psi_V=psi[jif][o].tolist())
    ysi = np.asarray(sol0['y'][-len(geo.y_si):], float); wsi = np.asarray(sol0['w_si'], float)[-len(geo.y_si):]
    k = np.clip(np.searchsorted(ysi, ys), 1, len(ysi) - 1)
    k = np.where(np.abs(ysi[k - 1] - ys) < np.abs(ysi[k] - ys), k - 1, k)
    assert np.max(np.abs(ysi[k] - ys)) < 1e-12, 'depth-line mapping'
    wn = wsi[k]
    xcol = np.array(sorted(set(np.round(xs, 16)))); ci = np.searchsorted(xcol, np.round(xs, 16))
    out['column_sheet_density'] = dict(x_cm=xcol.tolist(), N_cm2=np.bincount(ci, weights=n * wn, minlength=len(xcol)).tolist())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--deck', choices=['W', 'R'], required=True)
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--vg', type=float, required=True)
    ap.add_argument('--tox', type=float, default=20.0); ap.add_argument('--na', type=float, default=1e16)
    ap.add_argument('--hw', type=float, default=0.5); ap.add_argument('--lw', type=float, default=0.3)
    ap.add_argument('--ov', type=float, default=0.0, help='gate overlap (nm), deck R')
    ap.add_argument('--sigma', type=float, default=0.0, help='lateral erfc grading (nm), deck R')
    ap.add_argument('--psibi-shift', type=float, default=0.0, help='V, deck W')
    ap.add_argument('--devsim', action='store_true', help='also DEVSIM coupled-Jacobian cross-check')
    ap.add_argument('--tag', default='')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    name = (f"e7_{a.deck}_L{a.L:g}_lev{a.level}_VG{a.vg:g}_tox{a.tox:g}_NA{a.na:.0e}_hw{a.hw:g}"
            f"_ov{a.ov:g}_sg{a.sigma:g}_dpsi{a.psibi_shift:g}{('_' + a.tag) if a.tag else ''}.json")
    path = os.path.join(a.out, name)
    rec = dict(version=VERSION, start_utc=now(), args=vars(a), N_D=ND, lombardi=LOMB)
    t0 = time.time()
    set_bulk(a.na)
    rec['bulk'] = dict(N_A=D.NA, P0=D.P0, psi_b0=D.PSI_B0, W_dep_cm=MM.wdep_max_cm(na_cm3=a.na))
    geo = Geo(a)
    rec['geometry'] = geo.describe()
    rec['mesh'] = build_mesh(geo)
    D.build_physics()
    sol0, rec['equilibrium'] = equilibrium(geo, a)
    rec['t_equilibrium_s'] = time.time() - t0
    rec['result'] = analyse(geo, sol0)
    if a.devsim:
        dv = WP.devsim_G(); rec['devsim'] = dv
        rec['cross_check'] = (rec['result']['const_mu']['G_net'] - dv['G_E']) / dv['G_E']
    rec['wall_s'] = time.time() - t0
    rec['end_utc'] = now()
    json.dump(rec, open(path, 'w'), indent=1)
    c, m = rec['result']['const_mu'], rec['result']['lombardi_mu']
    print(f"dL_net {c['dL_net_nm']:.4f} slice {c['dL_slice_nm']:.4f} meshref {c['dL_net_meshref_nm']:.4f} "
          f"centre_vs_1col {c['centre_vs_1col']:.2e} | mu: dL_net {m['dL_net_nm']:.3f} slice {m['dL_slice_nm']:.3f} "
          f"mu[{m['mu_min']:.0f},{m['mu_max']:.0f}] | bounds {c['bounds_ok']},{m['bounds_ok']} | wall {rec['wall_s']:.0f}s "
          f"| cross {rec.get('cross_check')}")


if __name__ == '__main__':
    main()
