"""R-I(b) boundary inspection (WS4): read the actual equations at the boundary nodes from the
exported Newton system instead of assuming them.

For an R-I(b) device (lateral psi Neumann, phi_n Dirichlet), V_D = 0, given V_G, report for
selected Si/oxide nodes: psi - psi_1COL, psi - psi_b0 (bottom), row nnz, whether the bulk row
was dropped/moved by a contact or interface (exported permutation), and the columns coupled.
Nodes: Si corners (0,0), (L,0); bottom corners (0,xmax), (L,xmax); bottom neighbour (dx,xmax);
oxide corners (0,0), (L,0). Also: where in the domain the lateral non-uniformity sits.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import devsim as ds

from phase3 import runlog
from phase3.c1 import deck as D
from phase3.c1.c1b_jacobian import export
from phase3.c1.corner_study import eqnums
from phase3.c3.t1 import Geo, equilibrate, onecol_sol, _lut


def describe(eq, A, prow, region, var, nd):
    rows = eq[(region, var)]; r = int(rows[nd]); Ar = A.tocsr()
    cols = Ar.indices[Ar.indptr[r]:Ar.indptr[r + 1]]
    names = []
    for c in cols:
        for k, v in eq.items():
            hit = np.where(v == c)[0]
            if len(hit):
                names.append(f'{k[0]}:{k[1]}#{int(hit[0])}')
    return dict(row=r, nnz=int(len(cols)), bulk_row_permuted_to=int(prow[r]), cols=names[:12])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--L', type=float, default=1.0)
    ap.add_argument('--level', type=int, default=0); ap.add_argument('--vg', type=float, default=2.0)
    ap.add_argument('--body-corners', type=int, default=0); ap.add_argument('--out', required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    geo = Geo(a.L, a.level)
    rec = runlog.header(vars(a))
    rec['mesh'] = D.build_mesh_explicit(geo.xs, geo.y_ox, geo.y_si, corner='tie', body_corners=bool(a.body_corners))
    D.build_physics()
    sol0, rec['equilibrium'] = equilibrate(geo, a.vg, 'b')
    lut = _lut(sol0)
    eq = eqnums()
    A, b, Ad, bd, prow, pcopy, keys = export()
    xs, ys = map(np.array, D.coords(D.SI)); xo, yo = map(np.array, D.coords(D.OX))
    psi = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Potential'))
    pso = np.array(ds.get_node_model_values(device=D.DEV, region=D.OX, name='Potential'))
    L, xm, dx = geo.L_cm, float(ys.max()), float(geo.xs[1] - geo.xs[0])
    sel = {'si_corner_src(0,0)': (0, 0), 'si_corner_drn(L,0)': (L, 0), 'bottom_corner_src(0,xmax)': (0, xm),
           'bottom_corner_drn(L,xmax)': (L, xm), 'bottom_nb(dx,xmax)': (dx, xm), 'bottom_mid(L/2,xmax)': (float(geo.xs[len(geo.xs) // 2]), xm)}
    out = {}
    for name, (x, y) in sel.items():
        nd = int(np.argmin((xs - x) ** 2 + (ys - y) ** 2))
        out[name] = dict(x=float(xs[nd]), y=float(ys[nd]), psi_minus_1col=float(psi[nd] - lut[round(float(ys[nd]), 15)]),
                         psi_minus_psib0=float(psi[nd] - D.PSI_B0),
                         Potential=describe(eq, A, prow, D.SI, 'Potential', nd),
                         Electrons=describe(eq, A, prow, D.SI, 'Electrons', nd))
    for name, x in (('ox_corner_src(0,0)', 0.0), ('ox_corner_drn(L,0)', L)):
        nd = int(np.argmin((xo - x) ** 2 + yo ** 2))
        out[name] = dict(psi_minus_1col=float(pso[nd] - lut[0.0]), Potential=describe(eq, A, prow, D.OX, 'Potential', nd))
    # where is the lateral non-uniformity: deviation from 1COL, maximum per row, top 6 nodes
    dev = np.array([psi[i] - lut[round(float(ys[i]), 15)] for i in range(len(psi))])
    top = np.argsort(-np.abs(dev))[:8]
    rec['largest_deviation_from_1col'] = [dict(x=float(xs[i]), y=float(ys[i]), dev_V=float(dev[i])) for i in top]
    rec['nodes'] = out
    rec['contacts'] = {c: ds.get_contact_equation_list(device=D.DEV, contact=c) for c in ds.get_contact_list(device=D.DEV)}
    runlog.write_json(os.path.join(a.out, f'rib_inspect_L{a.L:g}_lev{a.level}_VG{a.vg}_bc{a.body_corners}.json'), rec)
    print('INSPECT', json.dumps(dict(top=rec['largest_deviation_from_1col'][:4],
                                     bottom={k: (v['psi_minus_psib0'], v['Potential']['nnz'], v['Potential']['bulk_row_permuted_to'])
                                             for k, v in out.items() if 'bottom' in k},
                                     contacts=rec['contacts'])))
    return 0


if __name__ == '__main__':
    sys.exit(main())
