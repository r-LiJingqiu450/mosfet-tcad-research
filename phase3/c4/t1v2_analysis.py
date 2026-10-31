"""Aggregate the registered T1-v2 runs (C-4 Task B/C/D). Verdicts are read from the run files
(computed by the registered rule); nothing is re-evaluated here. Also compares the historical
decks with the v1 records (results/phase3/C3/t1): x_U must reproduce v1 G_L1 and the AC path
must reproduce v1 Re Y (same solver, same state).
Usage: python phase3/c4/t1v2_analysis.py --in DIR [DIR ...] --v1 results/phase3/C3/t1 --out OUT
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def rows_from(paths):
    rows = []
    for p in paths:
        r = json.load(open(p)); cfg = r['config']
        for key, pt in r['points'].items():
            ev = pt['eval']; d = pt['diagnostics']
            row = dict(model=cfg['model'], lres=cfg['lres'] if cfg['model'] == 'c' else 0.0, L_um=cfg['L_um'],
                       level=cfg['level'], VG_V=pt['vg'], n_unknowns=ev['n_unknowns'], file=os.path.basename(p),
                       G_GCA_h=pt['G_GCA_h_S_per_cm'], G_M0=pt['M0_G0_S_per_cm'])
            for X in ('L1v2', 'L2v2'):
                for k, v in pt[X].items():
                    row[f'{X}_{k}'] = v
            row.update(eta_E=ev['E']['eta'], eta_M=ev['M']['eta'], eta_AC=ev['AC_f2']['eta'], eta_NS=ev['NS']['eta'],
                       r0_M=ev['M']['r0'], r0_NS=ev['NS']['r0'], kappa_G=ev['AC_f2']['kappa_G'],
                       consistency_ok=d['consistency']['ok'], consistency_ratio=d['consistency']['diff'] / d['consistency']['limit'],
                       B3v2_ok=d['B3v2']['ok'], B5v2_lin=d['B5v2']['lin'], B5v2_ok=d['B5v2']['ok'],
                       NS_vs_M=d['NS_vs_M_rel'], AC_vs_M=d['AC_vs_M_rel'],
                       lam_rel_diff=ev['adjoint']['lam_rel_diff'])
            if 'U' in ev:
                row.update(r0_U=ev['U']['r0'], eta_U=ev['U']['eta'], dwr_pred_U=d['dwr_validation_U']['predicted_lam_M'],
                           dwr_actual_U=d['dwr_validation_U']['actual_U_minus_M'])
            if 'decoupling' in pt:
                row.update(decoupling_rel=pt['decoupling']['rel_change'], decoupling_ok=pt['decoupling']['ok'])
            row.update(continuity_jump=pt['continuity']['max_abs_jump_V'],
                       gauss=min(pt['gauss']['closure_plus'], pt['gauss']['closure_minus']),
                       single_entry_oxide_rows=len(pt['oxide_lateral_single_entry_rows']),
                       psi_spread=pt['uniformity'][0]['max_spread'], n_spread=pt['uniformity'][1]['max_spread'],
                       bottom_max_psi_dev=max(abs(b['psi_minus_psib0']) for b in pt['bottom_rows']),
                       bottom_corner_rows=';'.join(f"{b['column']}:{b['Potential_row_nnz']}/{b['Electrons_row_nnz']}"
                                                   for b in pt['bottom_rows']),
                       wall_s=pt['wall_s'])
            rows.append(row)
    return rows


def v1_compare(rows, v1dir):
    out = []
    for r in rows:
        if r['model'] not in ('a', 'b_c3'):
            continue
        f = os.path.join(v1dir, f"t1_L{r['L_um']:g}_lev{r['level']}.json")
        if not os.path.exists(f):
            continue
        v1 = json.load(open(f))['points'][f"{'a' if r['model'] == 'a' else 'b'}:VG={r['VG_V']}"]
        out.append(dict(model=r['model'], L_um=r['L_um'], level=r['level'], VG_V=r['VG_V'],
                        v1_L1_verdict=v1['L1']['verdict'], v1_L2_verdict=v1['L2']['verdict'],
                        v1_r0_L1=v1['L1']['r0'], v2_r0_U=r.get('r0_U'), U_reproduces_v1=(r.get('r0_U') == v1['L1']['r0']),
                        v1_r0_L2=v1['L2']['r0'], v2_r0_L2=r['L2v2_r0'], L2_reproduces_v1=(r['L2v2_r0'] == v1['L2']['r0']),
                        v1_E_L1=v1['L1']['E'], v1_E_L2=v1['L2']['E'], v2_E_L1=r['L1v2_E'], v2_E_L2=r['L2v2_E'],
                        v2_L1_verdict=r['L1v2_verdict'], v2_L2_verdict=r['L2v2_verdict']))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--in', dest='inp', nargs='+', required=True)
    ap.add_argument('--v1', required=True); ap.add_argument('--out', required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    paths = sorted(p for d in a.inp for p in glob.glob(os.path.join(d, 't1v2_*.json')))
    rows = rows_from(paths)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(a.out, 't1v2_table.csv'), 'w', newline='') as fh:
        fh.write('# T1 v2 (registered t1v2-2026-10-30); r0 and E dimensionless; E are estimates, not bounds\n')
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    comp = v1_compare(rows, a.v1)
    if comp:
        with open(os.path.join(a.out, 'v1_vs_v2.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(comp[0].keys())); w.writeheader(); w.writerows(comp)
    summ = dict(n_points=len(rows), n_files=len(paths), by_model={})
    grp = collections.defaultdict(list)
    for r in rows:
        grp[(r['model'], r['lres'])].append(r)
    for (m, lr), R in sorted(grp.items()):
        s = dict(n=len(R))
        for X in ('L1v2', 'L2v2'):
            s[X] = dict(collections.Counter(r[f'{X}_verdict'] for r in R))
            s[X + '_wc'] = dict(collections.Counter(r[f'{X}_verdict_wc'] for r in R))
            s[X + '_max_abs_r0'] = max(abs(r[f'{X}_r0']) for r in R)
            s[X + '_max_E'] = max(r[f'{X}_E'] for r in R)
            s[X + '_max_E_wc'] = max(r[f'{X}_E_wc'] for r in R)
        s.update(consistency_fail=sum(not r['consistency_ok'] for r in R),
                 max_consistency_ratio=max(r['consistency_ratio'] for r in R),
                 B3v2_fail=sum(not r['B3v2_ok'] for r in R), max_eta_E=max(r['eta_E'] for r in R),
                 B5v2_fail=sum(not r['B5v2_ok'] for r in R), max_B5v2_lin=max(r['B5v2_lin'] for r in R),
                 max_abs_NS_vs_M=max(abs(r['NS_vs_M']) for r in R), max_abs_AC_vs_M=max(abs(r['AC_vs_M']) for r in R),
                 max_continuity=max(r['continuity_jump'] for r in R), max_gauss=max(r['gauss'] for r in R),
                 single_entry=sorted({r['single_entry_oxide_rows'] for r in R}),
                 max_psi_spread=max(r['psi_spread'] for r in R), max_bottom_psi_dev=max(r['bottom_max_psi_dev'] for r in R),
                 bottom_corner_rows=sorted({r['bottom_corner_rows'] for r in R if r['level'] == 0}),
                 max_kappa_G=max(r['kappa_G'] for r in R))
        dec = [r for r in R if 'decoupling_ok' in r]
        if dec:
            s['decoupling_fail'] = [dict(L_um=r['L_um'], level=r['level'], VG_V=r['VG_V'], rel=r['decoupling_rel'],
                                         E_L1v2=r['L1v2_E']) for r in dec if not r['decoupling_ok']]
            s['max_abs_decoupling'] = max(abs(r['decoupling_rel']) for r in dec)
        U = [r for r in R if 'r0_U' in r]
        if U:
            ratios = [r['dwr_pred_U'] / r['dwr_actual_U'] for r in U if r['dwr_actual_U'] != 0]
            s['dwr_validation_U'] = dict(n=len(ratios), min_ratio=min(ratios), max_ratio=max(ratios),
                                         median_ratio=sorted(ratios)[len(ratios) // 2],
                                         max_abs_actual=max(abs(r['dwr_actual_U']) for r in U))
        s['non_pass'] = [dict(L_um=r['L_um'], level=r['level'], VG_V=r['VG_V'],
                              L1v2=(r['L1v2_verdict'], r['L1v2_r0'], r['L1v2_E']), L2v2=(r['L2v2_verdict'], r['L2v2_r0'], r['L2v2_E']),
                              wc=(r['L1v2_verdict_wc'], r['L2v2_verdict_wc']))
                         for r in R if r['L1v2_verdict'] != 'PASS' or r['L2v2_verdict'] != 'PASS'
                         or r['L1v2_verdict_wc'] != 'PASS' or r['L2v2_verdict_wc'] != 'PASS']
        summ['by_model'][f'{m}' + (f'_lres{lr:g}' if m == 'c' else '')] = s
    if comp:
        summ['v1_reproduction'] = dict(n=len(comp), U_reproduces_v1=sum(c['U_reproduces_v1'] for c in comp),
                                       L2_reproduces_v1=sum(c['L2_reproduces_v1'] for c in comp),
                                       v1_non_pass=[c for c in comp if c['v1_L1_verdict'] != 'PASS' or c['v1_L2_verdict'] != 'PASS'])
    json.dump(summ, open(os.path.join(a.out, 't1v2_summary.json'), 'w'), indent=1, default=str)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
        cols = {'a': 'C0', 'b_c3': 'C3', 'b': 'C1', 'c': 'C2'}
        for k, X in enumerate(('L2v2', 'L1v2')):
            for m in cols:
                R = [r for r in rows if r['model'] == m]
                if not R:
                    continue
                ax[k].scatter([r['n_unknowns'] for r in R], [max(abs(r[f'{X}_r0']), 1e-16) for r in R], s=14,
                              color=cols[m], label=f'|r0| {m}')
                ax[k].scatter([r['n_unknowns'] for r in R], [3 * r[f'{X}_E'] for r in R], marker='_', s=60, color=cols[m])
            ax[k].axhline(1e-6, ls='--', c='k', lw=.8)
            ax[k].set(xscale='log', yscale='log', xlabel='unknowns', title=f'T1 v2 {X}: |r0| (dots) and 3E (bars)')
        ax[0].set_ylabel('relative'); ax[0].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, 'T1v2_overview.png'), dpi=150); fig.savefig(os.path.join(a.out, 'T1v2_overview.svg'))
    except Exception as exc:
        summ['figure_error'] = repr(exc); json.dump(summ, open(os.path.join(a.out, 't1v2_summary.json'), 'w'), indent=1, default=str)
    print(json.dumps({k: v for k, v in summ.items() if k != 'by_model'}, default=str)[:2000])
    for k, s in summ['by_model'].items():
        print(k, json.dumps({kk: s[kk] for kk in ('n', 'L1v2', 'L2v2', 'L1v2_wc', 'L2v2_wc', 'L1v2_max_E', 'L2v2_max_E',
                                                  'consistency_fail', 'B3v2_fail', 'B5v2_fail', 'max_gauss', 'single_entry',
                                                  'max_psi_spread', 'max_bottom_psi_dev')}, default=str))
    return 0


if __name__ == '__main__':
    sys.exit(main())
