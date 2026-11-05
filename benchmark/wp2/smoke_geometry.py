#!/usr/bin/env python3
"""WP2 smoke test (pre-registration section 5): with uniform conductivity, the network solvers of P1 (DEVSIM mesh
geometry + rdeck_common.terminal_solve) and P4 (fv_rdeck.Grid + fv_rdeck.terminal_solve) must reproduce the planning
factors of benchmark/specs/support/body_geometry.py to 1e-3 relative at matched grid spacing (h = 5 nm).

Geometry exactly as body_geometry.py: p region of unit conductivity; n+ wells (optionally enlarged by W_j) are
perfect conductors (fixed nodes) at drain = 1, source = 0; bottom y = 1.2 um fixed at 0 (B1) or free (B0).
Reported: G_DD, G_SD, G_DB in units of sigma_bulk for both paths, the planning values, and an FV-consistent
variant of the planning stencil (boundary edges with half weight, as box integration gives) for diagnosis.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

H_UM = 0.005
CASES = (('dev1', 1.0, 0.0), ('dev1', 1.0, 0.037), ('dev2', 2.0, 0.0), ('dev2', 2.0, 0.351))


def planning_solve():
    """Load solve() from the frozen planning script without running its module-level loop."""
    src = open(os.path.join(ROOT, 'benchmark/specs/support/body_geometry.py')).read()
    ns = {}
    exec(compile(src.split('\nfor name, L, wj in')[0], 'body_geometry.py', 'exec'), ns)
    return ns['solve']


def planning_fv(L, lw, hw, depth, wj, ohmic_body, h=0.005):
    """Planning stencil with box-integration boundary weights (edges on the outer boundary: weight 1/2)."""
    x = np.arange(-lw, L + lw + h / 2, h); y = np.arange(0, depth + h / 2, h)
    nx, ny = len(x), len(y); X, Y = np.meshgrid(x, y, indexing='ij')
    src = (X <= 0 + wj + 1e-12) & (Y <= hw + wj + 1e-12)
    drn = (X >= L - wj - 1e-12) & (Y <= hw + wj + 1e-12)
    body = (np.abs(Y - depth) < 1e-12) & ~src & ~drn if ohmic_body else np.zeros_like(src)
    idx = np.arange(nx * ny).reshape(nx, ny)
    wh = np.ones((nx - 1, ny)); wh[:, 0] = 0.5; wh[:, -1] = 0.5
    wv = np.ones((nx, ny - 1)); wv[0, :] = 0.5; wv[-1, :] = 0.5
    r = np.r_[idx[:-1, :].ravel(), idx[:, :-1].ravel()]; c = np.r_[idx[1:, :].ravel(), idx[:, 1:].ravel()]
    w = np.r_[wh.ravel(), wv.ravel()]
    W = sp.coo_matrix((np.r_[w, w], (np.r_[r, c], np.r_[c, r])), shape=(nx * ny,) * 2).tocsr()
    Lap = (sp.diags(np.asarray(W.sum(1)).ravel()) - W).tocsr()
    fixed = (src | drn | body).ravel(); phi = np.zeros(nx * ny); phi[drn.ravel()] = 1.0
    fr = np.where(~fixed)[0]; fx = np.where(fixed)[0]
    phi[fr] = spla.spsolve(Lap[fr][:, fr].tocsc(), -Lap[fr][:, fx] @ phi[fx])
    I = Lap @ phi
    return dict(I_D=float(I[drn.ravel()].sum()), I_S=float(I[src.ravel()].sum()), I_B=float(I[body.ravel()].sum()))


def masks(X, Y, L, wj, depth=1.2):
    src = (X <= 0 + wj + 1e-12) & (Y <= 0.5 + wj + 1e-12)
    drn = (X >= L - wj - 1e-12) & (Y <= 0.5 + wj + 1e-12)
    body = (np.abs(Y - depth) < 1e-12) & ~src & ~drn
    return src, drn, body


def p1_path(L, wj):
    import rdeck_common as RC
    ds = RC.ds
    h = H_UM * 1e-4
    xs = np.arange(-0.3e-4, (L + 0.3) * 1e-4 + h / 2, h); ys = np.arange(0, 1.2e-4 + h / 2, h)
    nx, ny = len(xs), len(ys)
    nid = lambda i, j: j * nx + i
    coords = []
    for j in range(ny):
        for i in range(nx):
            coords += [float(xs[i]), float(ys[j]), 0.0]
    el = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            el += [2, 0, a, b, c, 2, 0, a, c, d]
    mesh, dev = f'smoke_{L}_{wj}', f'smoke_dev_{L}_{wj}'
    ds.create_gmsh_mesh(mesh=mesh, coordinates=coords, elements=el, physical_names=['si'])
    ds.add_gmsh_region(gmsh_name='si', mesh=mesh, region='si', material='Si')
    ds.finalize_mesh(mesh=mesh)
    ds.create_device(mesh=mesh, device=dev)
    ds.edge_from_node_model(device=dev, region='si', node_model='node_index')
    n0i = np.array(ds.get_edge_model_values(device=dev, region='si', name='node_index@n0')).astype(int)
    n1i = np.array(ds.get_edge_model_values(device=dev, region='si', name='node_index@n1')).astype(int)
    c = np.array(ds.get_edge_model_values(device=dev, region='si', name='EdgeCouple')) * \
        np.array(ds.get_edge_model_values(device=dev, region='si', name='EdgeInverseLength'))
    X = np.array(ds.get_node_model_values(device=dev, region='si', name='x')) * 1e4
    Y = np.array(ds.get_node_model_values(device=dev, region='si', name='y')) * 1e4
    src, drn, body = masks(X, Y, L, wj)
    Lap = RC.laplacian(len(X), n0i, n1i, c)
    out = {}
    for b in ('B0', 'B1'):
        T = dict(source=np.where(src)[0], drain=np.where(drn)[0])
        if b == 'B1':
            T['body'] = np.where(body)[0]
        r = RC.terminal_solve(Lap, n0i, n1i, c, T, 'drain')
        out[b] = dict(G_DD=r['energy'], G_SD=-r['I']['source'], G_DB=-r['I'].get('body', 0.0))
    return out


def p4_path(L, wj):
    import fv_rdeck as F
    h = H_UM * 1e-4
    xs = np.arange(-0.3e-4, (L + 0.3) * 1e-4 + h / 2, h); ys = np.arange(0, 1.2e-4 + h / 2, h)
    g = F.Grid(xs, np.array([0.0]), ys)
    m = g.ksi > 0
    a, b, c = g.ea[m], g.eb[m], g.ksi[m]            # uniform conductivity 1
    src, drn, body = masks(g.X * 1e4, g.Y * 1e4, L, wj)
    out = {}
    for bd in ('B0', 'B1'):
        T = dict(source=np.where(src)[0], drain=np.where(drn)[0])
        if bd == 'B1':
            T['body'] = np.where(body)[0]
        r = F.terminal_solve(g.N, a, b, c, T, 'drain')
        out[bd] = dict(G_DD=r['energy'], G_SD=-r['I']['source'], G_DB=-r['I'].get('body', 0.0))
    return out


def main():
    solve = planning_solve()
    rows, ok_plan, ok_fv, ok_p1p4 = [], True, True, True
    for dev, L, wj in CASES:
        plan = {'B1': solve(L, 0.3, 0.5, 1.2, wj, True), 'B0': solve(L, 0.3, 0.5, 1.2, wj, False)}
        pfv = {'B1': planning_fv(L, 0.3, 0.5, 1.2, wj, True), 'B0': planning_fv(L, 0.3, 0.5, 1.2, wj, False)}
        p1, p4 = p1_path(L, wj), p4_path(L, wj)
        for b in ('B0', 'B1'):
            ref = dict(G_DD=plan[b]['I_D'], G_SD=-plan[b]['I_S'], G_DB=-plan[b]['I_B'])
            reffv = dict(G_DD=pfv[b]['I_D'], G_SD=-pfv[b]['I_S'], G_DB=-pfv[b]['I_B'])
            for q in ('G_DD', 'G_SD') + (('G_DB',) if b == 'B1' else ()):
                e1, e4 = p1[b][q] / ref[q] - 1, p4[b][q] / ref[q] - 1
                f1, f4 = p1[b][q] / reffv[q] - 1, p4[b][q] / reffv[q] - 1
                x14 = p1[b][q] / p4[b][q] - 1
                ok_plan &= abs(e1) <= 1e-3 and abs(e4) <= 1e-3
                ok_fv &= abs(f1) <= 1e-3 and abs(f4) <= 1e-3
                ok_p1p4 &= abs(x14) <= 1e-9
                rows.append(dict(device=dev, W_j_um=wj, body=b, quantity=q, planning=ref[q], planning_fv=reffv[q],
                                 P1=p1[b][q], P4=p4[b][q], P1_vs_planning=e1, P4_vs_planning=e4,
                                 P1_vs_planning_fv=f1, P4_vs_planning_fv=f4, P1_vs_P4=x14))
                print(f"{dev} wj={wj:5.3f} {b} {q}: plan={ref[q]:.6f} planFV={reffv[q]:.6f} P1={p1[b][q]:.6f} P4={p4[b][q]:.6f} "
                      f"| P1/plan-1={e1:+.2e} P4/plan-1={e4:+.2e} | P1/planFV-1={f1:+.2e} P4/planFV-1={f4:+.2e} | P1/P4-1={x14:+.1e}", flush=True)
    res = dict(criterion='|P/planning - 1| <= 1e-3 for P1 and P4 (pre-registration section 5)',
               pass_vs_planning_script=bool(ok_plan), pass_vs_planning_fv_weights=bool(ok_fv),
               P1_equals_P4_within_1e9=bool(ok_p1p4), rows=rows)
    out = os.path.join(ROOT, 'benchmark/wp2/results/smoke_geometry.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(res, f, indent=1)
    print('PASS vs planning script:', ok_plan, '| PASS vs FV-weight planning stencil:', ok_fv, '| P1 == P4:', ok_p1p4)
    return 0


if __name__ == '__main__':
    sys.exit(main())
