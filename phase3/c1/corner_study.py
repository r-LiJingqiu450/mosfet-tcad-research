"""Corner/interface (triple-point) defect study: original vs alternative constructions.

One configuration per process (DEVSIM state is global). Per V_G in {0.8, 1.25, 2.0} V
at V_D = 0 (R-I(a)-type lateral data from the independent 1COL solve unless stated):

  identity PREVIEW r0 = G_2D/G_GCA,h - 1 via L1 (exported matrix) and L2 (AC 0.1 Hz)
      -- a diagnostic, NOT a T1 result (no error budget, no D-2/D-3 analysis);
  decoupling: change of G_L1 when the lateral psi data follow V_D (s_couple = 1);
  potential continuity at every coincident oxide/Si node pair (incl. corners);
  which oxide nodes carry a Dirichlet (contact) equation besides the gate;
  oxide-corner equation row in the exported Newton matrix;
  Gauss closure: sum of contact charges vs total space charge;
  lateral uniformity of the V_D = 0 solution.

Configurations:
  A0  built-in mesher, lateral data from the 2-D Neumann run   (C-1 attempt 1)
  A0b built-in mesher, 1COL lateral data, no oxide treatment   (C-1 attempt 2)
  A1  built-in mesher, 1COL data, oxide node ADJACENT to the triple point pinned
      (C-1 attempt 3, 'legacy pin')
  X0  explicit mesh, no corner treatment (control)
  X1  explicit mesh, one-node oxide contact ON the triple point = continuity tie
  X2  explicit mesh, Si lateral contacts exclude the interface node (endpoint separation)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import scipy.sparse.linalg as spla

import devsim as ds

from phase3 import meshmap as MM, runlog
from phase3.c1 import deck as D, measure as M
from phase3.c1.c1b_jacobian import REC, recording_callback, export, snapshot, restore

VGS = (0.8, 1.25, 2.0)


def explicit_coords(level, n_lat0=20):
    xs = MM.lateral_map(1.0, n_lat0).nodes(level)
    y_si = MM.f0_depth_map().nodes(level)
    y_ox = MM.oxide_map(20.0, 5).nodes(level)
    return xs, y_ox, y_si


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


def g_l1(eq):
    A, b0, *_ = export()
    lu = spla.splu(A.tocsc())
    d = 1e-3
    ds.circuit_alter(name='V1', value=d)
    _, bd, *_ = export()
    ds.circuit_alter(name='V1', value=0.)
    s = lu.solve(bd - b0)
    return -s[int(eq[('circuit', 'V1.I')][0])] / d, A


def g_l2(f=0.1):
    ds.solve(type='ac', frequency=f, solver_type='direct')
    return -ds.get_circuit_node_value(solution='ssac_real', node='V1.I')


def continuity():
    xs, ys = map(np.array, D.coords(D.SI)); xo, yo = map(np.array, D.coords(D.OX))
    ps = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Potential'))
    po = np.array(ds.get_node_model_values(device=D.DEV, region=D.OX, name='Potential'))
    lat = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='PsiLat'))
    si_if = {round(float(x), 14): i for i, (x, y) in enumerate(zip(xs, ys)) if y == 0.0}
    out, worst = [], 0.0
    for j, (x, y) in enumerate(zip(xo, yo)):
        if y == 0.0 and round(float(x), 14) in si_if:
            i = si_if[round(float(x), 14)]
            d = float(po[j] - ps[i]); worst = max(worst, abs(d))
            if x == xo.min() or x == xo.max():
                out.append(dict(x_cm=float(x), psi_ox=float(po[j]), psi_si=float(ps[i]),
                                psi_lat_1col=float(lat[i]), jump_V=d))
    return dict(max_abs_jump_V=worst, corners=out)


def oxide_corner_row(eq, A):
    xo, yo = map(np.array, D.coords(D.OX))
    rows = eq[(D.OX, 'Potential')]
    Ar = A.tocsr()
    res = []
    for xc in (xo.min(), xo.max()):
        nd = int(np.where((xo == xc) & (yo == 0.0))[0][0])
        r = int(rows[nd]); cols = Ar.indices[Ar.indptr[r]:Ar.indptr[r + 1]]
        vals = Ar.data[Ar.indptr[r]:Ar.indptr[r + 1]]
        names = []
        for c in cols:
            for k, v in eq.items():
                hit = np.where(v == c)[0]
                if len(hit):
                    names.append(f'{k[0]}:{k[1]}#{int(hit[0])}')
        res.append(dict(x_cm=float(xc), row=r, nnz=int(len(cols)), cols=names,
                        vals=[float(v) for v in vals]))
    return res


def gauss_closure(eq):
    """Discrete Gauss law. Summing all bulk Poisson rows that are NOT replaced by a
    contact equation cancels interior fluxes; what remains is the displacement flux
    into the contact nodes plus the space charge of the non-contact nodes.
    Contact rows are identified from the exported permutation (row moved/dropped).
    Closure = |sum_c Q_c + s * sum_{non-contact Si} PNC*NodeVolume| / scale, s = +-1
    (DEVSIM's sign convention for get_contact_charge is established by the result)."""
    A, b, Ad, bd, prow, pcopy, keys = export()
    q = {}
    for c in ds.get_contact_list(device=D.DEV):
        try:
            q[c] = ds.get_contact_charge(device=D.DEV, contact=c, equation='PotentialEquation')
        except ds.error:
            pass
    rows = eq[(D.SI, 'Potential')]
    n_dev = int(min(v.min() for k, v in eq.items() if k[0] == 'circuit'))   # first circuit row
    # contact rows: bulk row dropped (-1, Dirichlet contact) or moved into a circuit row.
    # Rows moved into oxide rows are INTERFACE rows (fluxes combined), not contacts.
    is_contact = np.array([(prow[r] < 0) or (prow[r] >= n_dev) for r in rows])
    pnc = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='PotentialNodeCharge'))
    nv = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='NodeVolume'))
    rho_nc = float(np.sum((pnc * nv)[~is_contact])); rho_all = float(np.sum(pnc * nv))
    sq = float(sum(q.values()))
    scale = max(abs(v) for v in q.values())
    return dict(contact_charges_C_per_cm=q, sum_contact=sq, space_charge_noncontact=rho_nc,
                space_charge_all=rho_all, n_si_contact_rows=int(is_contact.sum()),
                closure_plus=abs(sq + rho_nc) / scale, closure_minus=abs(sq - rho_nc) / scale)


def oxide_lateral_rows(eq, A):
    """Rows of oxide nodes on the lateral edges (excluding the gate node): a node with
    an artificial Dirichlet would show a single-entry row."""
    xo, yo = map(np.array, D.coords(D.OX)); rows = eq[(D.OX, 'Potential')]; Ar = A.tocsr()
    out = []
    for xc in (xo.min(), xo.max()):
        idx = np.where(xo == xc)[0]; idx = idx[np.argsort(yo[idx])]
        for nd in idx[1:]:                          # skip gate node (top)
            r = int(rows[nd]); nnz = int(Ar.indptr[r + 1] - Ar.indptr[r])
            out.append(dict(x_cm=float(xc), y_cm=float(yo[nd]), row_nnz=nnz))
    single = [o for o in out if o['row_nnz'] == 1]
    return dict(n_lateral_oxide_nodes=len(out), single_entry_rows=single)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True, choices=['A0', 'A0b', 'A1', 'X0', 'X1', 'X2'])
    ap.add_argument('--level', type=int, default=0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    rec = runlog.header(dict(config=a.config, level=a.level, vgs=VGS, pm=D.pm_checklist()))
    t0 = time.time()
    if a.config in ('A0', 'A0b', 'A1'):
        rec['mesh'] = D.build_mesh(corner_contacts=(a.config == 'A1'))
    else:
        xs, y_ox, y_si = explicit_coords(a.level)
        rec['mesh'] = D.build_mesh_explicit(xs, y_ox, y_si,
                                            corner={'X0': 'none', 'X1': 'tie', 'X2': 'exclude'}[a.config])
    D.build_physics()
    rec['contacts'] = list(ds.get_contact_list(device=D.DEV))
    L = rec['mesh']['L_cm']
    rec['points'] = {}
    for vg in VGS:
        r = {}
        ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
        r['equilibrium'] = D.equilibrium(vg) if a.config == 'A0' else D.equilibrium_1col(vg)
        G = M.gca_conductance(L)['G_S_per_cm']
        eq = eqnums()
        gl1, A = g_l1(eq)
        gl2 = g_l2()
        r.update(G_GCA_h=G, G_L1=gl1, G_L2_0p1Hz=gl2, r0_L1=gl1 / G - 1, r0_L2=gl2 / G - 1,
                 L1_vs_L2=gl1 / gl2 - 1, continuity=continuity(),
                 uniformity=[M.lateral_uniformity('Potential'), M.lateral_uniformity('Electrons')],
                 oxide_corner_rows=oxide_corner_row(eq, A), gauss=gauss_closure(eq),
                 oxide_lateral=oxide_lateral_rows(eq, A),
                 n_unknowns=int(A.shape[0]))
        ds.set_parameter(device=D.DEV, name='s_couple', value=1.)
        D.converge()
        g1c, _ = g_l1(eqnums())
        r['decoupling_rel_change_L1'] = g1c / gl1 - 1
        ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
        D.converge()
        rec['points'][f'VG={vg}'] = r
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(out, f'{a.config}_L{a.level}.json'), rec)
    print('CORNER_DONE', a.config, a.level,
          json.dumps({k: [v['r0_L1'], v['continuity']['max_abs_jump_V']] for k, v in rec['points'].items()}))
    return 0


if __name__ == '__main__':
    sys.exit(main())
