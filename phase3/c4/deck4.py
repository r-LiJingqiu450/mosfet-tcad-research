"""C-4 reference decks: R-I(a) (unchanged), R-I(b) with corrected bottom corners, R-I(c)
reservoir strips. Builds on phase3/c1/deck.py (same physics, same contact equations for
gate/body/source/drain); only the mesh/contact geometry and two extra contacts differ.

Model geometry (DEVSIM coordinates: x along the channel, y = depth):

  'a'      R-I(a), C-3 deck, unchanged: lateral contacts source/drain on the full Si edge
           (rows interface..bottom) with psi = PsiLat (+ s_couple*bias) and phi_n Dirichlet;
           the bottom corners belong to the lateral contacts (their psi data = psi_b0 there).
  'b_c3'   R-I(b), C-3 deck, unchanged (historical): same mesh as 'a'; lateral contacts carry
           only phi_n. Open item R-I(b)-1: the bottom corner nodes carry no psi condition.
  'b'      R-I(b), corrected (C-4 Task C): lateral contacts on rows interface..bottom-1;
           each bottom corner is a separate one-node contact ('bs'/'bd') carrying BOTH
           psi = psi_b0 (bottom condition, R04 4.2) and phi_n = bias (lateral condition,
           R04 4.2 'over the full depth'). No node belongs to two contacts.
  'c'      R-I(c) reservoir strips (C-4 Task D; R04 4.3 'optional R-I(c) family'):
           Si strips of length l_res at each end, x in [-l_res, 0] and [L, L + l_res],
           same doping/physics; gate oxide and gate cover the strips. In each strip every
           node above the bottom row carries phi_n = bias (contact 'source'/'drain', electron
           equation only: psi natural = bulk Poisson row); the strip bottom row is contact
           'bs'/'bd' with psi = psi_b0 and phi_n = bias. Outer lateral boundaries (Si and
           oxide) are homogeneous Neumann. l_res -> 0 reproduces 'b' node for node.
           Channel: x in (0, L); the channel length L is the distance between the inner
           strip edges, as in 'a'/'b'.

All models keep the one-node oxide contacts oxc_s/oxc_d on the outer triple points; they carry
an equation only in 'a' (Research_08 section 5).
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import devsim as ds

from phase3 import meshmap as MM
from phase3.c1 import deck as D

N_RES0 = 10                       # strip cells at level 0 (uniform), N_k = round(10 r^k)
L_RES_FACTORS = (0.1, 0.3, 1.0)   # pre-specified (R04 4.3), times W_dep
STATE = {'model': None, 'dual': False}


def wdep_cm():
    return MM.wdep_max_cm()


def strip_nodes(l_res_cm, level):
    return MM.AxisMap('reservoir', (MM.Segment(0.0, l_res_cm, N_RES0, 'uniform'),)).nodes(level)


class Geo4:
    """Lateral/depth/oxide lines for one model, length and level."""

    def __init__(self, model, L_um, level, l_res_factor=0.0, n_lat0=50):
        self.model, self.L_um, self.level = model, L_um, level
        self.xs_ch = MM.lateral_map(L_um, n_lat0).nodes(level)
        self.L_cm = float(self.xs_ch[-1])
        self.y_si = MM.f0_depth_map().nodes(level)
        self.y_ox = MM.oxide_map(20.0, 5).nodes(level)
        self.l_res_factor = l_res_factor if model == 'c' else 0.0
        self.l_res_cm = self.l_res_factor * wdep_cm()
        if model == 'c':
            s = strip_nodes(self.l_res_cm, level)                       # 0..l_res
            self.xs = np.concatenate([s[:-1] - self.l_res_cm, self.xs_ch, self.L_cm + s[1:]])
            self.n_res = len(s) - 1
        else:
            self.xs = self.xs_ch.copy(); self.n_res = 0
        self.ic0 = self.n_res                       # column index of x = 0
        self.ic1 = self.ic0 + len(self.xs_ch) - 1   # column index of x = L
        assert self.xs[self.ic0] == 0.0 and self.xs[self.ic1] == self.L_cm

    def describe(self):
        return dict(model=self.model, L_um=self.L_um, level=self.level, l_res_factor=self.l_res_factor,
                    l_res_cm=self.l_res_cm, W_dep_cm=wdep_cm(), n_res_cells=self.n_res,
                    n_lateral_nodes=int(len(self.xs)), n_channel_cells=int(len(self.xs_ch) - 1),
                    x_min_cm=float(self.xs[0]), x_max_cm=float(self.xs[-1]),
                    strip_h_cm=(self.l_res_cm / self.n_res if self.n_res else None),
                    channel_h_cm=float(self.xs_ch[1] - self.xs_ch[0]))


def build_mesh(geo: Geo4):
    """Explicit tensor mesh (two right triangles per rectangle, as deck.build_mesh_explicit)."""
    model = geo.model
    xs = np.asarray(geo.xs, float); y_ox = np.asarray(geo.y_ox, float); y_si = np.asarray(geo.y_si, float)
    ys = np.concatenate([y_ox, y_si[1:]]); kif = len(y_ox) - 1
    nx, ny = len(xs), len(ys)
    nid = lambda i, j: j * nx + i
    coords = []
    for j in range(ny):
        for i in range(nx):
            coords += [float(xs[i]), float(ys[j]), 0.0]
    dual = model in ('b', 'c')
    names = [D.OX, D.SI, 'ox_si', 'gate', 'source', 'drain', 'body', 'oxc_s', 'oxc_d'] + (['bs', 'bd'] if dual else [])
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
    ic0, ic1 = geo.ic0, geo.ic1
    jlast = ny - 1 if not dual else ny - 2          # last row of the lateral/strip contacts
    for c, cols in (('source', range(0, ic0 + 1)), ('drain', range(ic1, nx))):
        for i in cols:
            for j in range(kif, jlast):
                el += [1, P[c], nid(i, j), nid(i, j + 1)]
    # body: channel bottom strictly inside (0, L)
    for i in range(ic0 + 1, ic1 - 1):
        el += [1, P['body'], nid(i, ny - 1), nid(i + 1, ny - 1)]
    if dual:
        for c, cols in (('bs', list(range(0, ic0 + 1))), ('bd', list(range(ic1, nx)))):
            if len(cols) == 1:
                el += [1, P[c], nid(cols[0], ny - 1), nid(cols[0], ny - 1)]
            for i0, i1 in zip(cols[:-1], cols[1:]):
                el += [1, P[c], nid(i0, ny - 1), nid(i1, ny - 1)]
    el += [1, P['oxc_s'], nid(0, kif), nid(0, kif)]
    el += [1, P['oxc_d'], nid(nx - 1, kif), nid(nx - 1, kif)]
    m = 'm4'
    ds.create_gmsh_mesh(mesh=m, coordinates=coords, elements=el, physical_names=names)
    ds.add_gmsh_region(gmsh_name=D.OX, mesh=m, region=D.OX, material='Ox')
    ds.add_gmsh_region(gmsh_name=D.SI, mesh=m, region=D.SI, material='Si')
    ds.add_gmsh_interface(gmsh_name='ox_si', mesh=m, name='ox_si', region0=D.OX, region1=D.SI)
    contacts = [('gate', D.OX), ('source', D.SI), ('drain', D.SI), ('body', D.SI), ('oxc_s', D.OX), ('oxc_d', D.OX)]
    if dual:
        contacts += [('bs', D.SI), ('bd', D.SI)]
    for c, r in contacts:
        ds.add_gmsh_contact(gmsh_name=c, mesh=m, name=c, region=r, material='metal')
    ds.finalize_mesh(mesh=m)
    ds.create_device(mesh=m, device=D.DEV)
    D.CORNER['on'] = True; D.CORNER['mode'] = 'tie'
    STATE.update(model=model, dual=dual)
    return dict(kind='explicit_c4', nx=nx, ny=ny, kif=kif, ic0=ic0, ic1=ic1,
                lateral_contact_rows=[kif, jlast], dual_bottom_contacts=dual, **geo.describe())


def add_dual_contacts():
    """'bs'/'bd': psi = psi_b0 (no circuit node, like the body) and phi_n = bias
    (drain side on circuit node drain_bias, so its electron current is drain current)."""
    for c, b in (('bs', 'source_bias'), ('bd', 'drain_bias')):
        ds.contact_node_model(device=D.DEV, contact=c, name=f'{c}_psi', equation='Potential - psi_b0')
        ds.contact_node_model(device=D.DEV, contact=c, name=f'{c}_psi:Potential', equation='1')
        ds.contact_equation(device=D.DEV, contact=c, name='PotentialEquation', node_model=f'{c}_psi',
                            edge_charge_model='contactcharge_edge')
        eq = f'Electrons - n_i*exp((Potential - {b})/V_t)'
        ds.contact_node_model(device=D.DEV, contact=c, name=f'{c}_n', equation=eq)
        ds.contact_node_model(device=D.DEV, contact=c, name=f'{c}_n:Electrons', equation='1')
        ds.contact_node_model(device=D.DEV, contact=c, name=f'{c}_n:Potential',
                              equation=f'-n_i*exp((Potential - {b})/V_t)/V_t')
        kw = dict(device=D.DEV, contact=c, name='ElectronContinuityEquation', node_model=f'{c}_n',
                  edge_current_model='ElectronCurrent')
        if c == 'bd':
            ds.contact_node_model(device=D.DEV, contact=c, name=f'{c}_n:{b}',
                                  equation=f'n_i*exp((Potential - {b})/V_t)/V_t')
            kw['circuit_node'] = b
        ds.contact_equation(**kw)


def build(geo: Geo4):
    if geo.model in ('a', 'b_c3'):
        # historical C-3 deck, byte-for-byte the same builder call as phase3/c3/t1.py
        rec = D.build_mesh_explicit(geo.xs, geo.y_ox, geo.y_si, corner='tie')
        rec.update(kind='explicit_c3', dual_bottom_contacts=False, **geo.describe())
        STATE.update(model=geo.model, dual=False)
        D.build_physics()
        return rec
    rec = build_mesh(geo)
    D.build_physics()
    if STATE['dual']:
        add_dual_contacts()
    return rec


def mode_of(model):
    """psi treatment of the lateral contacts: Dirichlet only in R-I(a)."""
    return 'a' if model == 'a' else 'b'
