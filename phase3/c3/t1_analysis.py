"""Aggregate the C-3/C-4 T1 runs: verdict table (pre-registered rule, unchanged), error
components, boundary diagnostics, mesh convergence of G_h vs M0 (R04 7.4 with ACTUAL h_rep
ratios), L3 cross-check and small-delta behaviour of r(delta). Writes CSV, JSON, figures.
Usage: python phase3/c3/t1_analysis.py --in results/phase3/C3/t1 --out results/phase3/C3/analysis
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np


def observed_order(d_fine, d_coarse, r_fine, r_coarse, it=200):
    """Observed order from three solutions with NON-constant refinement ratios.
    d_fine = f1 - f2, d_coarse = f2 - f3 (same sign, monotone), r_fine = h2/h1, r_coarse = h3/h2.
    For f_k = c + a h_k^p:  d_coarse/d_fine = r_fine^p (r_coarse^p - 1)/(r_fine^p - 1), so
      p = [ln(d_coarse/d_fine) + ln((r_fine^p - 1)/(r_coarse^p - 1))] / ln r_fine
    (Celik et al., J. Fluids Eng. 130, 078001, 2008, monotone case), solved by fixed-point
    iteration. With r_fine == r_coarse it reduces to ln(d_coarse/d_fine)/ln r.
    v1 of this analysis (commit 896d9d9) used only the first term, which biases p by
    ln((r_fine^p-1)/(r_coarse^p-1))/ln r_fine (about -0.008 for the D-1 depth ratios)."""
    lr = math.log(d_coarse / d_fine)
    p = lr / math.log(r_fine)
    for _ in range(it):
        pn = (lr + math.log((r_fine ** p - 1) / (r_coarse ** p - 1))) / math.log(r_fine)
        if abs(pn - p) < 1e-14:
            return pn
        p = pn
    return p


def classify(f, h):
    """R04 7.4 with actual ratios. f, h: finest first (f1 finest)."""
    f1, f2, f3, f4 = f
    r12, r23, r34 = h[1] / h[0], h[2] / h[1], h[3] / h[2]
    d12, d23, d34 = f1 - f2, f2 - f3, f3 - f4
    out = dict(R123=None, R234=None, p123=None, p234=None, cls=None, U=None)
    if d12 == 0 and d23 == 0:
        out.update(cls='converged_to_roundoff', U=0.0); return out
    R123 = d12 / d23 if d23 != 0 else math.inf
    R234 = d23 / d34 if d34 != 0 else math.inf
    out.update(R123=R123, R234=R234)
    if 0 < R123 < 1 and 0 < R234 < 1:
        p123 = observed_order(d12, d23, r12, r23); p234 = observed_order(d23, d34, r23, r34)
        out.update(p123=p123, p234=p234)
        if 0.5 <= p123 <= 3 and abs(p123 - p234) <= 0.25 * p123:
            out.update(cls='asymptotic', U=1.25 * abs(d12) / (r12 ** p123 - 1),
                       richardson=f1 + d12 / (r12 ** p123 - 1)); return out
        out.update(cls='monotone_nonasymptotic', U=3 * max(abs(x - f1) for x in f)); return out
    out.update(cls='oscillatory', U=3 * (max(f[:3]) - min(f[:3])) / 2); return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--out', required=True); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    runs = {}
    for p in sorted(glob.glob(os.path.join(a.inp, 't1_L*_lev*.json'))):
        r = json.load(open(p)); runs[(r['config']['L_um'], r['config']['level'])] = r
    rows = []
    for (L, lev), r in sorted(runs.items()):
        for key, p in r['points'].items():
            mode, vg = p['mode'], p['vg']
            row = dict(L_um=L, level=lev, mode=mode, VG_V=vg, n_unknowns=p['n_unknowns'],
                       G_GCA_h_S_per_um=p['G_GCA_h_S_per_cm'] * 1e-4, G_M0_S_per_um=p['M0_G0_S_per_cm'] * 1e-4,
                       eq_newton_its=p['equilibrium']['newton_iterations'])
            for X in ('L1', 'L2'):
                for c in ('r0', 'E_freq', 'E_eq', 'E_lin', 'E', 'verdict'):
                    if c in p[X]:
                        row[f'{X}_{c}'] = p[X][c]
            row['decoupling_rel'] = p.get('decoupling', {}).get('rel_change')
            row['decoupling_ok'] = p.get('decoupling', {}).get('ok')
            l3 = p.get('L3')
            if l3:
                row.update(L3_status=l3['status'], L3_r0=l3.get('r0'), L3_E=l3.get('E'),
                           L3_consistent=(l3.get('consistency_with_L2') or {}).get('ok'),
                           r_at_8mV=l3['points'][-1]['r'], c2_r_over_delta2=l3['points'][0]['r'] / l3['points'][0]['delta_V'] ** 2)
            row.update(continuity_jump_V=p['continuity']['max_abs_jump_V'],
                       gauss=min(p['gauss']['closure_plus'], p['gauss']['closure_minus']),
                       single_entry_oxide_rows=len(p['oxide_lateral_single_entry_rows']),
                       lateral_psi_spread_V=p['uniformity'][0]['max_spread'],
                       lateral_n_spread_rel=p['uniformity'][1]['max_spread'],
                       wall_s=p['wall_s'])
            rows.append(row)
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0].keys()).index(k) if k in rows[0] else 999)
    with open(os.path.join(a.out, 't1_table.csv'), 'w', newline='') as fh:
        fh.write('# T1 table; G in S/um; r0, E dimensionless (estimates, not bounds)\n')
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    # verdict summary
    summ = {}
    for X in ('L1', 'L2'):
        v = [r[f'{X}_verdict'] for r in rows]
        summ[X] = {k: v.count(k) for k in set(v)}
        summ[X + '_max_abs_r0'] = max(abs(r[f'{X}_r0']) for r in rows)
        summ[X + '_max_E'] = max(r[f'{X}_E'] for r in rows)
    summ['n_points'] = len(rows)
    summ['max_continuity_jump_V'] = max(r['continuity_jump_V'] for r in rows)
    summ['max_gauss'] = max(r['gauss'] for r in rows)
    summ['single_entry_rows_by_mode'] = {m: sorted({r['single_entry_oxide_rows'] for r in rows if r['mode'] == m}) for m in 'ab'}
    summ['max_lateral_psi_spread_V_by_mode'] = {m: max(r['lateral_psi_spread_V'] for r in rows if r['mode'] == m) for m in 'ab'}
    summ['decoupling_all_ok'] = all(r['decoupling_ok'] for r in rows if r['mode'] == 'a')
    summ['decoupling_failures'] = [dict(L_um=r['L_um'], level=r['level'], VG_V=r['VG_V'],
                                        rel_change=r['decoupling_rel'], E_L1=r['L1_E'],
                                        ratio=abs(r['decoupling_rel']) / r['L1_E'])
                                   for r in rows if r['mode'] == 'a' and not r['decoupling_ok']]
    summ['non_pass_points'] = [dict(L_um=r['L_um'], level=r['level'], mode=r['mode'], VG_V=r['VG_V'],
                                    L1=(r['L1_verdict'], r['L1_r0'], r['L1_E']), L2=(r['L2_verdict'], r['L2_r0'], r['L2_E']))
                               for r in rows if r['L1_verdict'] != 'PASS' or r['L2_verdict'] != 'PASS']
    summ['max_decoupling_rel'] = max(abs(r['decoupling_rel']) for r in rows if r['mode'] == 'a')
    l3rows = [r for r in rows if r.get('L3_status')]
    summ['L3'] = dict(n=len(l3rows), statuses=sorted({r['L3_status'] for r in l3rows}),
                      all_consistent=all(r['L3_consistent'] for r in l3rows if r['L3_status'] == 'ok'))
    # mesh convergence of G_h (per L, VG; identical for both modes since G_GCA,h is the 1COL value)
    conv = []
    from phase3 import meshmap as MM
    dm = MM.f0_depth_map(); lat = {}
    for (L, vg) in sorted({(r['L_um'], r['VG_V']) for r in rows}):
        levs = sorted({r['level'] for r in rows if r['L_um'] == L and r['VG_V'] == vg})
        if levs != [0, 1, 2, 3]:
            continue
        G = {r['level']: r['G_GCA_h_S_per_um'] for r in rows if r['L_um'] == L and r['VG_V'] == vg and r['mode'] == 'a'}
        GL2 = {r['level']: (1 + r['L2_r0']) * r['G_GCA_h_S_per_um'] for r in rows
               if r['L_um'] == L and r['VG_V'] == vg and r['mode'] == 'a'}
        h = [1.2e-4 / (len(dm.nodes(k)) - 1) for k in (3, 2, 1, 0)]
        G0 = next(r['G_M0_S_per_um'] for r in rows if r['L_um'] == L and r['VG_V'] == vg)
        c = classify([G[3], G[2], G[1], G[0]], h)
        c2 = classify([GL2[3], GL2[2], GL2[1], GL2[0]], h)
        conv.append(dict(L_um=L, VG_V=vg, G_M0=G0, G_levels=[G[k] for k in range(4)], rel_diff_finest=G[3] / G0 - 1,
                         cls=c['cls'], p123=c['p123'], p234=c['p234'], U_disc_rel=(c['U'] or 0) / G[3],
                         T2form_pass=abs(G[3] / G0 - 1) <= 2 * (c['U'] or 0) / G[3] if c['cls'] == 'asymptotic' else None,
                         richardson_rel_minus_M0=(c.get('richardson', math.nan) / G0 - 1),
                         G_L2_cls=c2['cls'], G_L2_p123=c2['p123']))
    with open(os.path.join(a.out, 'G_mesh_convergence.csv'), 'w', newline='') as fh:
        fh.write('# G = channel conductance at V_D = 0; S/um; U_disc estimate per R04 7.4 with actual h_rep ratios\n')
        if conv:
            w = csv.DictWriter(fh, fieldnames=list(conv[0].keys())); w.writeheader(); w.writerows(conv)
    summ['G_convergence'] = conv
    json.dump(summ, open(os.path.join(a.out, 't1_summary.json'), 'w'), indent=1, default=str)
    # figures
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for X, mk in (('L1', 'o'), ('L2', 's')):
            for m, col in (('a', 'C0'), ('b', 'C3')):
                sel = [r for r in rows if r['mode'] == m]
                ax[0].scatter([r['n_unknowns'] for r in sel], [abs(r[f'{X}_r0']) for r in sel], marker=mk,
                              facecolors='none' if X == 'L1' else col, edgecolors=col, label=f'|r0| {X} R-I({m})')
                ax[0].scatter([r['n_unknowns'] for r in sel], [3 * r[f'{X}_E'] for r in sel], marker='_', color=col, s=80)
        ax[0].axhline(1e-6, ls='--', c='k', lw=.8); ax[0].text(2e4, 1.3e-6, 'E_req = 1e-6', fontsize=8)
        ax[0].set(xscale='log', yscale='log', xlabel='unknowns', ylabel='|r0| (markers), 3E (bars)', title='T1: identity residual vs 3E')
        ax[0].legend(fontsize=7, ncol=2)
        l3rows_ = [r for r in rows if r.get('L3_status') == 'ok']
        for m, col in (('a', 'C0'), ('b', 'C3')):
            for L in sorted({r['L_um'] for r in l3rows_}):
                sel = sorted([r for r in l3rows_ if r['mode'] == m and r['L_um'] == L], key=lambda r: r['VG_V'])
                ax[1].plot([r['VG_V'] for r in sel], [abs(r['c2_r_over_delta2']) for r in sel], marker='o' if m == 'a' else 's',
                           color=col, alpha=.4 + .2 * math.log10(L + 1), label=f'R-I({m}) L={L:g} um')
        ax[1].set(yscale='log', xlabel='V_G (V)', ylabel='|r(delta)|/delta^2 at 0.25 mV (V^-2)',
                  title='Finite-delta ratio r(delta) ~ c delta^2 (level 0)')
        ax[1].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, 'T1_overview.png'), dpi=150); fig.savefig(os.path.join(a.out, 'T1_overview.svg'))
    except Exception as exc:
        summ['figure_error'] = repr(exc)
        json.dump(summ, open(os.path.join(a.out, 't1_summary.json'), 'w'), indent=1, default=str)
    print(json.dumps({k: v for k, v in summ.items() if k != 'G_convergence'}, indent=1, default=str))
    for c in conv:
        print('CONV', c['L_um'], c['VG_V'], c['cls'], c['p123'], 'rel %.3e U %.3e T2form %s rich %.2e' % (c['rel_diff_finest'], c['U_disc_rel'], c['T2form_pass'], c['richardson_rel_minus_M0']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
