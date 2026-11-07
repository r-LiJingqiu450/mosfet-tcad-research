#!/usr/bin/env python3
"""WP8, grid convergence of the A-R references (pre-registration benchmark/specs/WP8_*.md, section 4).

For each A-R state and each of G_DD and G_SD, the P1 network values on the three published meshes of the device
(levels 0, 1, 2; read from the released case files) are analysed with the three-grid procedure of Celik et al. (2008),
which generalizes Roache's grid convergence index (GCI) to non-constant refinement ratios:

  h_i ~ N_i^(-1/2) (N_i: number of mesh nodes, nx * ny), r21 = h2/h1, r32 = h3/h2 (1 = finest = level 2);
  eps21 = f2 - f1, eps32 = f3 - f2, s = sign(eps32/eps21);
  p = |ln|eps32/eps21| + q(p)| / ln r21, q(p) = ln((r21^p - s)/(r32^p - s)), by fixed-point iteration from q = 0;
  f_ext = (r21^p f1 - f2)/(r21^p - 1), e_a21 = |(f1 - f2)/f1|, GCI_fine = 1.25 e_a21/(r21^p - 1);
  asymptotic ratio GCI_32/(r21^p GCI_21) (close to 1 in the asymptotic range).

Also reported: the registered tolerance of the case, and the independent finite-volume value (P4, its refined mesh)
against f_ext. Writes benchmark/wp8/results/gci_wp8.json; prints a Markdown table.

Usage (repository root):  python benchmark/wp8/gci_wp8.py
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
B = os.path.join(ROOT, 'benchmark')
OUT = os.path.join(HERE, 'results')
FS = 1.25


def three_grid(f1, f2, f3, r21, r32, max_it=200):
    """Celik et al. (2008) procedure. f1 finest. Returns a dict (p None if it cannot be determined)."""
    e21, e32 = f2 - f1, f3 - f2
    out = dict(f1=f1, f2=f2, f3=f3, r21=r21, r32=r32, eps21=e21, eps32=e32)
    if e21 == 0 or e32 == 0:
        out.update(s=0, monotone=False, p=None, note='zero difference')
        return out
    s = 1.0 if e32 / e21 > 0 else -1.0
    q, p = 0.0, None
    for _ in range(max_it):
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        arg = (r21 ** p_new - s) / (r32 ** p_new - s)
        if arg <= 0:
            p = None
            break
        q_new = math.log(arg)
        if p is not None and abs(p_new - p) < 1e-12:
            p = p_new
            break
        p, q = p_new, q_new
    out.update(s=int(s), monotone=s > 0, p=p)
    if p is None or p <= 0:
        return out
    rp = r21 ** p
    f_ext = (rp * f1 - f2) / (rp - 1)
    e_a21 = abs((f1 - f2) / f1)
    gci21 = FS * e_a21 / (rp - 1)
    e_a32 = abs((f2 - f3) / f2)
    gci32 = FS * e_a32 / (r32 ** p - 1)
    out.update(f_ext=f_ext, rel_f1_vs_ext=f1 / f_ext - 1, e_a21=e_a21, e_ext21=abs((f_ext - f1) / f_ext), GCI_fine=gci21,
               GCI_32=gci32, asymptotic_ratio=gci32 / (rp * gci21))
    return out


def mesh_nodes(mesh_id):
    m = json.load(open(os.path.join(B, 'meshes', f'mesh_{mesh_id}.json')))
    return m['nx'] * m['ny']


def tolerance(case, q):
    """The registered relative tolerance of the case for q (uncertainty.tolerance.<q>_rel)."""
    return case['uncertainty']['tolerance'][f'{q}_rel']


def analyse():
    rows = []
    for p in sorted(glob.glob(os.path.join(B, 'cases', 'AR__*.json'))):
        c = json.load(open(p))
        smr = c['observables']['same_mesh_references']
        by_level = {v['level']: (mid, v) for mid, v in smr.items()}
        if sorted(by_level) != [0, 1, 2]:
            continue
        N = {lev: mesh_nodes(by_level[lev][0]) for lev in (0, 1, 2)}
        r21, r32 = math.sqrt(N[2] / N[1]), math.sqrt(N[1] / N[0])
        p4 = {x['name']: x for x in c['cross_paths']}.get('independent_finite_volume_refined', {})
        for q, key in (('G_DD', 'G_DD_S_per_cm'), ('G_SD', 'G_SD_S_per_cm')):
            f1, f2, f3 = (by_level[lev][1][key] for lev in (2, 1, 0))
            g = three_grid(f1, f2, f3, r21, r32)
            tol = tolerance(c, q)
            g.update(state=c['device']['state_id'], case=c['id'], quantity=q, tolerance=tol, nodes=N,
                     ext_within_tolerance=(abs(g['rel_f1_vs_ext']) <= tol) if g.get('f_ext') is not None else None,
                     P4=p4.get(key), P4_rel_vs_ext=(p4[key] / g['f_ext'] - 1) if p4.get(key) and g.get('f_ext') else None)
            rows.append(g)
    return rows


def main():
    rows = analyse()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'gci_wp8.json'), 'w') as f:
        json.dump(dict(procedure='Celik et al. (2008) three-grid GCI, Fs = 1.25, h ~ N^-1/2', rows=rows), f, indent=1)
    print('| State | quantity | f(L0) | f(L1) | f(L2) | r21, r32 | monotone | p | f_ext | f(L2)/f_ext - 1 | GCI_fine | '
          'asympt. ratio | tolerance | within | P4/f_ext - 1 |')
    print('|---|---|---:|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---|---:|')
    for g in rows:
        def e(v, fmt='{:+.2e}'):
            return 'n/a' if v is None else fmt.format(v)
        print(f"| {g['state']} | {g['quantity']} | {g['f3']:.6e} | {g['f2']:.6e} | {g['f1']:.6e} | {g['r21']:.3f}, {g['r32']:.3f} | "
              f"{'yes' if g['monotone'] else 'no'} | {e(g.get('p'), '{:.2f}')} | {e(g.get('f_ext'), '{:.6e}')} | "
              f"{e(g.get('rel_f1_vs_ext'))} | {e(g.get('GCI_fine'), '{:.2e}')} | {e(g.get('asymptotic_ratio'), '{:.3f}')} | "
              f"{g['tolerance']:.2e} | {g['ext_within_tolerance']} | {e(g.get('P4_rel_vs_ext'))} |")
    return 0


if __name__ == '__main__':
    sys.exit(main())
