#!/usr/bin/env python3
"""WP8 exploratory analysis X: summary of the XS_*.json records (pre-registration section 8). Descriptive only.

Writes benchmark/wp8/results/exploratory/summary_xsens.json and prints Markdown tables.
Usage (repository root):  python benchmark/wp8/xsens_summary.py
"""
from __future__ import annotations

import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, 'results', 'exploratory')
PROD = os.path.join(HERE, 'results', 'production')


def load():
    out = []
    for p in sorted(glob.glob(os.path.join(EXP, 'XS_*.json'))):
        r = json.load(open(p))
        ss = json.load(open(os.path.join(PROD, os.path.basename(p).replace('XS_', 'SS_'))))
        out.append((r, ss))
    return out


def device(r):
    a = r['args']
    return f"L={a['L']:g} um, N_A={a['na']:.0e}, t_ox={a['tox']:g} nm"


def main():
    rows = []
    for r, ss in load():
        X2, X3 = r['X2'], r['X3']
        G_net = r['G_net']
        rows.append(dict(
            state=r['args']['id'], device=device(r), n=r['system']['n'], G_net=G_net,
            x_vs_lambda=r['solutions']['rel_x_vs_lambda'], x_vs_M1=r['solutions']['rel_x_vs_M1'],
            omega_x=r['solutions']['x_omega_history'][-1], omega_lambda=r['solutions']['lambda_omega_history'][-1],
            kappa_A=X2['kappa_A'], abs_rms_u=X2['abs_rms_u'], abs_sum_u=X2['abs_sum_u'], rel_rms_u=X2['rel_rms_u'],
            unrefined_abs_rms_u=X2['unrefined']['abs_rms_u'], unrefined_kappa_vs_G_net=X2['unrefined']['kappa_A_vs_G_net'],
            kappa_G_wp2=X2['unrefined']['kappa_G_wp2'],
            abs_err_IR=ss['M1']['G'] - G_net, e_IR=abs(ss['M1']['rel_vs_G_net']), S=ss['M2']['spread_rel_G_net'],
            shifts=[x['measured_shift'] for x in r['X1']], dG1=[x['dG1'] for x in r['X1']], dG2=[x['dG2'] for x in r['X1']],
            rel_dG1=[x['rel_dG1_vs_measured'] for x in r['X1']], rel_dG12=[x['rel_dG12_vs_measured'] for x in r['X1']],
            rho_rms=[x['rho_rms'] for x in r['X1']], ulp_mean=[x['ulp_mean'] for x in r['X1']],
            same_as_M2=[x.get('same_as_M2_record') for x in r['X1']],
            by_class=X3['by_class'], by_row_block=X3['by_row_block'], by_column_block=X3['by_column_block'],
            by_block_and_class=X3['by_block_and_class'], top=X3['top_rows'][:3], wall_s=r['wall_s']))
    with open(os.path.join(EXP, 'summary_xsens.json'), 'w') as f:
        json.dump(rows, f, indent=1)

    print('### X2: first-order sensitivity to the Jacobian entries (refined x and lambda)\n')
    print('| State | G_net (S/cm) | kappa_A | u (sum c^2)^1/2 (S/cm) | u sum|c| (S/cm) | same, unrefined | WP2 kappa_G | '
          '|G_IR - G_net| (S/cm) | x vs lambda |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for x in rows:
        print(f"| {x['state']} | {x['G_net']:.3e} | {x['kappa_A']:.2e} | {x['abs_rms_u']:.2e} | {x['abs_sum_u']:.2e} | "
              f"{x['unrefined_abs_rms_u']:.2e} | {x['kappa_G_wp2']:.1e} | {abs(x['abs_err_IR']):.2e} | {x['x_vs_lambda']:+.1e} |")
    print('\n### X1: response to the M2 perturbations (S/cm)\n')
    print('| State | seed | measured G_k - G_IR | first order dG1 | second order dG2 | dG1/measured - 1 | (dG1+dG2)/measured - 1 | '
          'rho rms | mean (ulp) | = M2 record |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|---|')
    for x in rows:
        for k in range(len(x['shifts'])):
            print(f"| {x['state']} | {k} | {x['shifts'][k]:+.4e} | {x['dG1'][k]:+.4e} | {x['dG2'][k]:+.2e} | {x['rel_dG1'][k]:+.2e} | "
                  f"{x['rel_dG12'][k]:+.2e} | {x['rho_rms'][k]:.3e} | {x['ulp_mean'][k]:+.4f} | {x['same_as_M2'][k]} |")
    print('\n### X3: share of sum|c_ij| by the class of the row node\n')
    classes = sorted({k for x in rows for k in x['by_class']})
    print('| State | ' + ' | '.join(classes) + ' |')
    print('|---|' + '---:|' * len(classes))
    for x in rows:
        print(f"| {x['state']} | " + ' | '.join(f"{x['by_class'].get(c, 0):.3f}" for c in classes) + ' |')
    print('\n### X3: share by equation block (rows / columns)\n')
    blocks = ('si_potential', 'si_electrons', 'ox_potential', 'circuit')
    print('| State | ' + ' | '.join(f'{b} rows' for b in blocks) + ' | ' + ' | '.join(f'{b} cols' for b in blocks) + ' |')
    print('|---|' + '---:|' * 8)
    for x in rows:
        print(f"| {x['state']} | " + ' | '.join(f"{x['by_row_block'].get(b, 0):.3f}" for b in blocks) + ' | '
              + ' | '.join(f"{x['by_column_block'].get(b, 0):.3f}" for b in blocks) + ' |')
    return 0


if __name__ == '__main__':
    sys.exit(main())
