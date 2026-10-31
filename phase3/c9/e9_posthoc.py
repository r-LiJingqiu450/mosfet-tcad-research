"""C-9 post-hoc analysis (NOT pre-registered; reported as such in the C-9 report, erratum E9-1 and addendum A9).

1. Mesh diagnostic for the halo devices: production level 0 vs level 1, and the short-device end mesh (a = 0.175 um)
   applied to L = 1 um (diag9_mesh.py), constant and Lombardi mobility, probe9 v1 (Lombardi impurity bug) and v2.
2. Affine-identity check of the maximum-g_m threshold shift: V_T(L) computed from the L = 1 um solves alone,
   G(L) = sigma_inf / (L - 2 dL(V_G)), compared with the direct 2-D solves at L = 0.5 and 2 um.
3. Level-1 strong-halo threshold shifts from the affine identity (mesh sensitivity of the reported dV_T).
Output: results/phase3/C9/analysis/e9_posthoc.json
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np  # noqa: E402

from phase3.c9 import e9_analysis as E  # noqa: E402

R9 = 'results/phase3/C9'
OUT = os.path.join(R9, 'analysis')
SRC = {'prod_v1': f'{R9}/runs', 'prod_v2': f'{R9}/runs_v2', 'lev1_v1': f'{R9}/diag/level1', 'lev1_v2': f'{R9}/diag/level1_v2',
       'aend_v1': f'{R9}/diag/aend0.175', 'aend_v2': f'{R9}/diag/aend0.175_v2'}


def load_dir(d):
    out = {}
    for p in glob.glob(os.path.join(d, 'e9_*.json')):
        r = json.load(open(p))
        if 'end_utc' not in r:
            continue
        a = r['args']; c, m = r['result']['const_mu'], r['result']['lombardi_mu']
        out[(a['deck'], a['halo_n'], a['L'], a['level'], a['vg'])] = dict(
            dL=c['dL_net_nm'], dLmu=m['dL_net_nm'], Gmu=m['G_net'], GGCAmu=m['G_GCA_1COL_mu'], version=r.get('version'))
    return out


def sweep(S, deck, hn, L, lev, f):
    ks = sorted(k for k in S if k[:4] == (deck, hn, L, lev))
    return np.array([k[4] for k in ks]), np.array([S[k][f] for k in ks])


def vt_from(v, G, vmax=2.0):
    return E.max_gm_vt(v, G, vmax)[0]


def affine_vt(S, deck, hn, lev, Ls, ref_L=1.0):
    v, GG = sweep(S, deck, hn, ref_L, lev, 'GGCAmu'); _, dl = sweep(S, deck, hn, ref_L, lev, 'dLmu')
    if len(v) < 5:
        return None
    vinf = vt_from(v, GG * ref_L)
    return {f'{L:g}': 1e3 * (vt_from(v, GG * ref_L / (L - 2e-3 * dl)) - vinf) for L in Ls}


def main():
    S = {k: load_dir(d) for k, d in SRC.items()}
    res = dict(sources=SRC, counts={k: len(v) for k, v in S.items()})
    # 1. mesh diagnostic at V_G = 1.32 V
    rows = []
    for tag, src, deck, hn, L, lev in (('production L=0.5 (a=0.175 um), level 0', 'prod', 'R', 2e18, 0.5, 0),
                                       ('production L=1 (a=0.35 um), level 0', 'prod', 'R', 2e18, 1.0, 0),
                                       ('production L=2 (a=0.35 um), level 0', 'prod', 'R', 2e18, 2.0, 0),
                                       ('L=1 on the short-device end mesh (a=0.175 um), level 0', 'aend', 'R', 2e18, 1.0, 0),
                                       ('L=0.5 (a=0.175 um), level 1', 'lev1', 'R', 2e18, 0.5, 1),
                                       ('L=1 (a=0.35 um), level 1', 'lev1', 'R', 2e18, 1.0, 1),
                                       ('halo 5e17, L=1, level 0', 'prod', 'R', 5e17, 1.0, 0),
                                       ('halo 5e17, L=1, level 1', 'lev1', 'R', 5e17, 1.0, 1),
                                       ('no halo, L=1, level 0', 'prod', 'R', 0.0, 1.0, 0),
                                       ('no halo, L=1, level 1', 'lev1', 'R', 0.0, 1.0, 1)):
        k = (deck, hn, L, lev, 1.32); row = dict(case=tag)
        for ver in ('v1', 'v2'):
            s = S.get(f'{src}_{ver}', {})
            if src == 'prod' and hn == 0.0:
                s = S['prod_v1']          # no-halo states are unaffected by erratum E9-1 (not re-run)
            if k in s:
                row[f'dL_const_{ver}'] = s[k]['dL']; row[f'dL_lombardi_{ver}'] = s[k]['dLmu']
        rows.append(row)
    res['mesh_diag_VG1.32'] = rows
    # level-1 vs level-0 sweeps at L = 1 um, strong halo
    sw = {}
    for ver in ('v1', 'v2'):
        p0 = S[f'prod_{ver}'] if ver == 'v2' else S['prod_v1']; p1 = S[f'lev1_{ver}']
        v1, a1 = sweep(p1, 'R', 2e18, 1.0, 1, 'dL'); _, m1 = sweep(p1, 'R', 2e18, 1.0, 1, 'dLmu')
        v0, a0 = sweep(p0, 'R', 2e18, 1.0, 0, 'dL'); _, m0 = sweep(p0, 'R', 2e18, 1.0, 0, 'dLmu')
        if len(v1) == len(v0) and len(v1) > 0 and np.allclose(v1, v0):
            sw[ver] = dict(VG=v1.tolist(), dL_const_lev0=a0.tolist(), dL_const_lev1=a1.tolist(), dL_lomb_lev0=m0.tolist(),
                           dL_lomb_lev1=m1.tolist(), max_rel_const=float(np.max(np.abs(a1 / a0 - 1))),
                           max_rel_lomb=float(np.max(np.abs(m1 / m0 - 1))))
    res['level_sweep_R_halo2e18_L1'] = sw
    # 2. affine identity for the threshold shift
    aff = {}
    for name, src, deck, hn in (('W_nohalo', 'prod_v1', 'W', 0.0), ('R_nohalo', 'prod_v1', 'R', 0.0),
                                ('R_halo5e17_v2', 'prod_v2', 'R', 5e17), ('R_halo2e18_v2', 'prod_v2', 'R', 2e18),
                                ('R_halo2e18_v1', 'prod_v1', 'R', 2e18)):
        s = S[src]; pred = affine_vt(s, deck, hn, 0, (0.5, 1.0, 2.0))
        if pred is None:
            continue
        direct = {}
        for L in (0.5, 1.0, 2.0):
            v, G = sweep(s, deck, hn, L, 0, 'Gmu'); _, GG = sweep(s, deck, hn, L, 0, 'GGCAmu')
            if len(v) >= 5:
                direct[f'{L:g}'] = 1e3 * (vt_from(v, G) - vt_from(v, GG))
        aff[name] = dict(dVT_direct_mV=direct, dVT_affine_from_L1_mV=pred)
    res['affine_identity_dVT'] = aff
    # 3. level-1 strong-halo threshold shifts (affine identity from the level-1 L = 1 um sweep)
    lv = {}
    for ver in ('v1', 'v2'):
        p = affine_vt(S[f'lev1_{ver}'], 'R', 2e18, 1, (0.5, 1.0, 2.0))
        if p is not None:
            lv[ver] = p
    res['dVT_halo2e18_level1_affine'] = lv
    json.dump(res, open(os.path.join(OUT, 'e9_posthoc.json'), 'w'), indent=1)
    print(json.dumps(res, indent=1)[:6000])


if __name__ == '__main__':
    main()
