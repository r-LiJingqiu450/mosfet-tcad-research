#!/usr/bin/env python3
"""Write the node lines of every WP2 deck-R mesh (A-R and held-out states) so that the meshes can be rebuilt
without this repository.

The meshes are the endprobe deck-R tensor-product maps used by P1/P2 (pre-registration section 4):
  lateral  phase3.c6.rprobe.lateral_nodes(L, level, l_w = 0.3 um)
  depth    phase3.c7.endprobe.depth_map(N_A).nodes(level)   (0 <= y <= 1.2 um)
  oxide    oxide_map(t_ox, 5).nodes(level)                  (-t_ox <= y <= 0)
Each rectangle (i, j)-(i+1, j+1) is split into the triangles [(i,j), (i+1,j), (i+1,j+1)] and [(i,j), (i+1,j+1), (i,j+1)]
(phase3/c7/endprobe.build_mesh). On these right triangles the box-integration coupling of the diagonal edge is zero,
so the Si network is the 5-point stencil on the node lines.

No device simulation is run. Each file is cross-checked against the mesh record of the matching P1 production run.

Usage (repository root):  python benchmark/wp2/mesh_lines.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rdeck_common as C  # noqa: E402

OUT = os.path.join(HERE, 'results', 'meshes')
PROD = os.path.join(HERE, 'results', 'production')

# (geometry key, N_A, t_ox nm, L um, levels, the states run on that mesh: every P1 record is cross-checked)
GEOMS = [('dev1_L1', 1e18, 2.0, 1.0, (0, 1, 2), ('AR-1', 'AR-2', 'AR-3', 'HO-1', 'HO-2')),
         ('dev2_L2', 1e16, 20.0, 2.0, (0, 1, 2), ('AR-4', 'AR-5', 'AR-6', 'HO-3', 'HO-4')),
         ('dev1_L2', 1e18, 2.0, 2.0, (1,), ('HO-5',))]


def lines(na, tox, L, level):
    a, geo = C.make_geo(L, level, tox, na, 0.0)
    xs = np.asarray(geo.xs, float)
    y_ox = np.asarray(geo.y_ox, float)
    y_si = np.asarray(geo.y_si, float)
    g0 = int(np.argmin(np.abs(xs + geo.ov_cm)))
    g1 = int(np.argmin(np.abs(xs - geo.L_cm - geo.ov_cm)))
    kif = len(y_ox) - 1
    return dict(
        units='cm', x_cm=[float(v) for v in xs], y_oxide_cm=[float(v) for v in y_ox], y_si_cm=[float(v) for v in y_si],
        nx=len(xs), ny=len(y_ox) + len(y_si) - 1, kif=kif, jw_si=int(geo.jw_si), contact_extent_cm=float(geo.h_w_cm),
        gate_columns=[g0, g1], x_junctions_cm=[float(xs[geo.i0]), float(xs[geo.i1])],
        rules=dict(
            node_lines='tensor product of x_cm with (y_oxide_cm followed by y_si_cm[1:]); y = 0 is the Si/SiO2 interface (shared line)',
            triangles='rectangle (i,j)-(i+1,j+1) split into [(i,j),(i+1,j),(i+1,j+1)] and [(i,j),(i+1,j+1),(i,j+1)]; '
                      'the diagonal edge has zero box-integration coupling, so the Si network is the 5-point stencil',
            source_contact='x = x_cm[0], silicon nodes with 0 <= y <= y_si_cm[jw_si] (= contact_extent_cm)',
            drain_contact='x = x_cm[-1], silicon nodes with 0 <= y <= y_si_cm[jw_si]',
            body_contact='y = y_si_cm[-1] = 1.2e-4 cm, all x',
            gate_contact='oxide top y = y_oxide_cm[0], columns gate_columns[0] .. gate_columns[1] (x = 0 .. L)',
            oxide_triple_points='the two Si/SiO2 nodes at x = x_cm[0] and x = x_cm[-1] are tied to the adjacent contact potential'),
        generator=dict(lateral=f'phase3.c6.rprobe.lateral_nodes(L={L:g}, level={level}, lw=0.3)',
                       depth=f'phase3.c7.endprobe.depth_map({na:.0e}).nodes({level})',
                       oxide=f'oxide_map({tox:g}, 5).nodes({level})'))


def main():
    os.makedirs(OUT, exist_ok=True)
    for key, na, tox, L, levels, sids in GEOMS:
        for lev in levels:
            d = lines(na, tox, L, lev)
            d.update(mesh_id=f'{key}_lev{lev}', N_A_cm3=na, t_ox_nm=tox, L_um=L, level=lev)
            blob = json.dumps(dict(x=d['x_cm'], yo=d['y_oxide_cm'], ys=d['y_si_cm']), sort_keys=True).encode()
            d['node_lines_sha256'] = hashlib.sha256(blob).hexdigest()
            chk = {}
            for sid, body in ((s_, b_) for s_ in sids for b_ in ('B1', 'B0')):
                p = os.path.join(PROD, f'P1_{sid}_{body}_lev{lev}.json')
                if os.path.exists(p):
                    r = json.load(open(p))
                    m, g = r['mesh'], r['geometry']
                    ok = (m['nx'] == d['nx'] and m['ny'] == d['ny'] and m['kif'] == d['kif'] and m['jw'] == d['kif'] + d['jw_si']
                          and list(m['gate_cols']) == d['gate_columns'] and g['contact_extent_cm'] == d['contact_extent_cm'])
                    chk[f'P1_{sid}_{body}_lev{lev}'] = bool(ok)
                    assert ok, (key, lev, body)
            d['cross_check_vs_P1_records'] = chk
            with open(os.path.join(OUT, f"mesh_{d['mesh_id']}.json"), 'w', encoding='utf-8', newline='\n') as f:
                f.write(json.dumps(d, sort_keys=True, indent=1) + '\n')
            print(d['mesh_id'], d['nx'], d['ny'], d['nx'] * d['ny'], 'nodes', chk)
    return 0


if __name__ == '__main__':
    sys.exit(main())
