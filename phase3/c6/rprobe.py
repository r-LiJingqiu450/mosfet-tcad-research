"""C-6 exploratory probe R1: zero-bias conductance with RESOLVED abrupt n+ wells (R-II type), to
test Ward's equipotential-well truncation used in wprobe.py. EXPLORATORY; F0 physics otherwise.

Geometry (DEVSIM x along the channel, y = depth):
  Si x in [-l_w, L + l_w], depth 0..1.2 um; abrupt box wells N_D = 1e20 for x < 0 or x > L and
  y <= h_w (nodes exactly at x = 0 or x = L are channel nodes, i.e. p-type); channel p-type N_A = 1e16.
  Oxide 20 nm over the whole width; gate contact on the oxide top over 0 <= x <= L only (no
  overlap: Ward's geometry); oxide top over the wells Neumann.
  Well contacts: far Si edges x = -l_w and x = L + l_w for 0 <= y <= h_w, ohmic:
  psi = V_T asinh((N_D - N_A)/(2 n_i)) and phi_n = 0 / V_D; oxide triple points tied to the same psi.
  Lateral Si edges below h_w: Neumann. Bottom psi = psi_b0, zero electron flux.
Equilibrium: seed = 1COL slice in the channel, charge-neutral psi in the wells; Newton (log-damped) with
a gate homotopy from V_G = 0 if the direct solve fails.
Conductance: same exact network / DEVSIM cross-check / bounds as wprobe.py.
dL_X = (L/2)(1 - G_GCA/G_X) with L the METALLURGICAL channel length; the well sheet resistance is
included in G (constant mu; it is O(1e-5) of the channel resistance here).
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

import devsim as ds

from phase3 import meshmap as MM
from phase3.c1 import deck as D
from phase3.c3.t1 import onecol_sol, g_gca, _lut
from phase3.c6 import wprobe as WP

ND = 1e20
PSI_W = D.VT * math.asinh((ND - D.NA) / (2 * D.NI))
VERSION = 'rprobe-2026-10-31'


def lateral_nodes(L_um, level, lw_um=0.3, a_um=0.35, h0_nm=0.5, hmax_nm=20.0):
    core = WP.lateral_map(L_um, a_um, h0_nm, hmax_nm).nodes(level)
    lw, h0, hm = lw_um * 1e-4, h0_nm * 1e-7, hmax_nm * 1e-7
    best = None
    for n in range(5, 4000):
        Q = MM.stretch_for(lw, h0, n)
        last = lw * (Q - Q ** ((n - 1) / n)) / (Q - 1)
        if last <= hm:
            best = (n, Q); break
    n, Q = best
    right = MM.Segment(0.0, lw, n, 'geometric', Q).nodes(level)          # fine at 0
    left = -right[::-1]
    L = float(core[-1])
    return np.concatenate([left[:-1], core, L + right[1:]]), len(left) - 1


class GeoR:
    def __init__(self, L_um, level, h_w_um=0.5, lw_um=0.3):
        self.L_um, self.level, self.lw_cm = L_um, level, lw_um * 1e-4
        self.xs, self.i0 = lateral_nodes(L_um, level, lw_um)
        self.i1 = len(self.xs) - 1 - self.i0
        self.L_cm = float(self.xs[self.i1] - self.xs[self.i0])
        assert self.xs[self.i0] == 0.0
        self.y_si = MM.f0_depth_map().nodes(level)
        self.y_ox = MM.oxide_map(20.0, 5).nodes(level)
        self.jw_si = int(np.argmin(np.abs(self.y_si - h_w_um * 1e-4)))
        self.h_w_cm = float(self.y_si[self.jw_si])

    def describe(self):
        return dict(L_um=self.L_um, level=self.level, l_w_cm=self.lw_cm, n_lateral=int(len(self.xs)),
                    n_depth_si=int(len(self.y_si)), h_w_cm=self.h_w_cm, i0=self.i0, i1=self.i1,
                    dx_min_nm=float(np.diff(self.xs).min() * 1e7), dx_max_nm=float(np.diff(self.xs).max() * 1e7))


def build_mesh(geo: GeoR):
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
    for i in range(geo.i0, geo.i1):
        el += [1, P['gate'], nid(i, 0), nid(i + 1, 0)]
    jw = kif + geo.jw_si
    for c, i in (('source', 0), ('drain', nx - 1)):
        for j in range(kif, jw):
            el += [1, P[c], nid(i, j), nid(i, j + 1)]
    for i in range(nx - 1):
        el += [1, P['body'], nid(i, ny - 1), nid(i + 1, ny - 1)]
    el += [1, P['oxc_s'], nid(0, kif), nid(0, kif)]
    el += [1, P['oxc_d'], nid(nx - 1, kif), nid(nx - 1, kif)]
    m = 'mr'
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


def well_mask(geo, xs, ys):
    return ((xs < -1e-15) | (xs > geo.L_cm * (1 + 1e-12))) & (ys <= geo.h_w_cm * (1 + 1e-12))


def set_doping(geo):
    xs, ys = map(np.array, D.coords(D.SI))
    nd = np.where(well_mask(geo, xs, ys), ND - D.NA, -D.NA)
    ds.node_solution(device=D.DEV, region=D.SI, name='NetDop')
    ds.set_node_values(device=D.DEV, region=D.SI, name='NetDop', values=list(nd))
    D._model(D.SI, 'PotentialNodeCharge', '-ElectronCharge*kahan3(Holes, -Electrons, NetDop)',
             ('Potential', 'Electrons'))
    return nd


def contact_mask(region, geo, xs, ys):
    x0, x1 = float(geo.xs[0]), float(geo.xs[-1])
    edge = (np.abs(xs - x0) < 1e-15) | (np.abs(xs - x1) < 1e-15)
    if region == D.SI:
        return edge & (ys <= geo.h_w_cm * (1 + 1e-12))
    return edge & (np.abs(ys) < 1e-15)


def equilibrium(geo, vg, tau=1e-12):
    sol0 = onecol_sol(geo, vg)
    lut = _lut(sol0)
    xs, ys = map(np.array, D.coords(D.SI)); xo, yo = map(np.array, D.coords(D.OX))
    wm = well_mask(geo, xs, ys)
    psi_si = np.array([lut[round(float(y), 15)] for y in ys]); psi_si[wm] = PSI_W
    psi_ox = np.array([lut[round(float(y), 15)] for y in yo])
    over_well = (xo < 0) | (xo > geo.L_cm)
    psi_ox[over_well] = PSI_W
    for r, v in ((D.SI, psi_si), (D.OX, psi_ox)):
        lat = v.copy()
        ds.set_node_values(device=D.DEV, region=r, name='Potential', values=list(v))
        pl = np.where(contact_mask(r, geo, *(map(np.array, D.coords(r)))), PSI_W, lat)
        ds.set_node_values(device=D.DEV, region=r, name='PsiLat', values=list(pl))
    ds.set_node_values(device=D.DEV, region=D.SI, name='Electrons', values=[D.NI * math.exp(p / D.VT) for p in psi_si])
    ds.set_parameter(device=D.DEV, name='gate_bias', value=vg)
    ds.circuit_alter(name='V1', value=0.)
    D.set_lateral_dirichlet(True)
    log = {}
    r = D.newton(tau, max_it=80, extra=2)
    log['direct'] = dict(converged=r['converged'], iterations=r['iterations'], last=r['history'][-3:])
    if not r['converged']:
        raise D.NotConverged(f'R1 equilibrium failed: {r["history"][-4:]}')
    return sol0, log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--vg', type=float, required=True); ap.add_argument('--hw', type=float, default=0.5)
    ap.add_argument('--lw', type=float, default=0.3)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f'r1_L{a.L:g}_lev{a.level}_VG{a.vg:g}_hw{a.hw:g}_lw{a.lw:g}.json')
    rec = dict(version=VERSION, start_utc=WP.now(), args=vars(a), psi_well_V=PSI_W, N_D=ND)
    t0 = time.time()
    geo = GeoR(a.L, a.level, a.hw, a.lw)
    rec['geometry'] = geo.describe()
    rec['mesh'] = build_mesh(geo)
    D.build_physics()
    nd = set_doping(geo)
    rec['n_well_nodes'] = int(np.sum(nd > 0))
    sol0, rec['equilibrium'] = equilibrium(geo, a.vg)
    # decomposition with the metallurgical L: reuse wprobe with a contact mask on the far edges
    WP.contact_mask = lambda region, g, xs, ys: contact_mask(region, geo, xs, ys)
    dc = WP.decomposition(geo, sol0)
    # wprobe's dl() used geo.L_cm (metallurgical) -> correct definition; surface profile keeps all x
    rec['decomposition'] = dc
    rec['devsim'] = WP.devsim_G()
    dv = rec['devsim']
    rec['cross_check'] = dict(net_vs_devsim_rel=(dc['G_net'] - dv['G_E']) / dv['G_E'],
                              devsim_E_vs_M_rel=(dv['G_E'] - dv['G_M']) / dv['G_E'])
    rec['wall_s'] = time.time() - t0
    rec['end_utc'] = WP.now()
    json.dump(rec, open(path, 'w'), indent=1)
    print({k: dc[k] for k in ('G_GCA', 'G_net', 'G_slice', 'G_rows', 'dL_net_nm', 'dL_slice_nm', 'dL_spread_nm', 'bounds_ok')})
    print('devsim', json.dumps(dv), 'cross', json.dumps(rec['cross_check']), 'wall', rec['wall_s'])


if __name__ == '__main__':
    main()
