#!/usr/bin/env python3
"""WP8 registered evaluation (pre-registration benchmark/specs/WP8_*.md, sections 3 and 4).

Small-signal mechanism, per state (records benchmark/wp8/results/production/SS_*.json):
  e_E, e_M      |G/G_net - 1| of the two double-precision solves (M0); e_dbl = min(e_E, e_M)
  e_IR          |G_IR/G_net - 1| of the exact solve of the same double-precision system (M1)
  S             max |G_k - G_IR| / |G_net| over the three one-ulp perturbations of the Jacobian (M2)
  P-M0  M0 reproduces WP2's G_E, G_M and G_net to <= 1e-12 relative (else the state is excluded)
  P-M1  among states with e_dbl > 1e-2: M1 converged and e_IR <= 1e-2 e_dbl in at least all but one
  P-M2  among the same states: S <= 1e-2 e_dbl in at least all but one
Grid convergence (benchmark/wp8/gci_wp8.py), 6 A-R states x {G_DD, G_SD}:
  P-G1  monotone convergence in >= 10 of 12 pairs
  P-G2  observed order 1 <= p <= 3 in every monotone pair
  P-G3  |f(L2)/f_ext - 1| <= the registered tolerance in all 12 pairs
Writes benchmark/wp8/results/evaluation_wp8.json; prints Markdown tables.

Usage (repository root):  python benchmark/wp8/evaluate_wp8.py
"""
from __future__ import annotations

import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gci_wp8 as GC  # noqa: E402

RES = os.path.join(HERE, 'results')
REPRO_TOL = 1e-12
SHARE = 1e-2
ERR_GATE = 1e-2


def mechanism():
    rows = []
    for p in sorted(glob.glob(os.path.join(RES, 'production', 'SS_*.json'))):
        r = json.load(open(p))
        m0, m1, m2 = r['M0'], r['M1'], r['M2']
        repro = all(abs(m0.get(k, 1.0)) <= REPRO_TOL for k in ('rel_vs_wp2_G_E', 'rel_vs_wp2_G_M')) and m0.get('same_G_net', False)
        e_E, e_M = abs(m0['rel_G_E']), abs(m0['rel_G_M'])
        e_dbl = min(e_E, e_M)
        e_IR = abs(m1['rel_vs_G_net'])
        S = r['M2']['spread_rel_G_net']
        m1b = r.get('M1b', {})
        rows.append(dict(state=r['args']['id'], n=r['system']['n'], G_net=r['G_net'], reproduces_wp2=repro,
                         rel_G_E=m0['rel_G_E'], rel_G_M=m0['rel_G_M'], e_dbl=e_dbl, M1_converged=m1['converged'],
                         M1_iterations=m1['iterations'], rel_G_IR=m1['rel_vs_G_net'], e_IR=e_IR,
                         omega_dd=m1['omega_dd_final'], S=S, M2_all_converged=m2['all_converged'],
                         M1b_rel_vs_G_IR=m1b.get('rel_vs_G_IR'), M1b_converged=m1b.get('converged'),
                         in_gate=e_dbl > ERR_GATE,
                         M1_ok=m1['converged'] and e_IR <= SHARE * e_dbl, M2_ok=S <= SHARE * e_dbl,
                         resolved=m1['converged'] and e_IR <= 1e-6, wall_s=r['wall_s']))
    valid = [x for x in rows if x['reproduces_wp2']]
    gate = [x for x in valid if x['in_gate']]
    n = len(gate)
    pm1 = sum(x['M1_ok'] for x in gate)
    pm2 = sum(x['M2_ok'] for x in gate)
    summary = dict(n_records=len(rows), n_reproducing=len(valid), excluded=[x['state'] for x in rows if not x['reproduces_wp2']],
                   n_in_gate=n, P_M1=dict(n_ok=pm1, n=n, supported=n > 0 and pm1 >= n - 1),
                   P_M2=dict(n_ok=pm2, n=n, supported=n > 0 and pm2 >= n - 1),
                   n_resolved_1e6=sum(x['resolved'] for x in valid))
    return rows, summary


def grid():
    rows = GC.analyse()
    mono = [g for g in rows if g['monotone']]
    p_ok = [g for g in mono if g.get('p') is not None and 1.0 <= g['p'] <= 3.0]
    within = [g for g in rows if g.get('ext_within_tolerance')]
    summary = dict(n_pairs=len(rows), P_G1=dict(n_monotone=len(mono), supported=len(mono) >= 10),
                   P_G2=dict(n_in_range=len(p_ok), n_monotone=len(mono), supported=len(p_ok) == len(mono) and len(mono) > 0),
                   P_G3=dict(n_within=len(within), n=len(rows), supported=len(within) == len(rows)))
    return rows, summary


def main():
    mrows, msum = mechanism()
    grows, gsum = grid()
    out = dict(mechanism=dict(rows=mrows, summary=msum), grid=dict(rows=grows, summary=gsum))
    with open(os.path.join(RES, 'evaluation_wp8.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print('### Small-signal mechanism (level 1, body B1)\n')
    print('| State | n | G_net (S/cm) | G_E rel | G_M rel | M1: G_IR rel | it | conv | M2 spread S | M1b vs M1 | gate | '
          'P-M1 | P-M2 |')
    print('|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|---|---|')
    for x in mrows:
        b = 'n/a' if x['M1b_rel_vs_G_IR'] is None else f"{x['M1b_rel_vs_G_IR']:+.1e}"
        print(f"| {x['state']} | {x['n']} | {x['G_net']:.4e} | {x['rel_G_E']:+.2e} | {x['rel_G_M']:+.2e} | {x['rel_G_IR']:+.2e} | "
              f"{x['M1_iterations']} | {x['M1_converged']} | {x['S']:.1e} | {b} | {'yes' if x['in_gate'] else 'no'} | "
              f"{x['M1_ok'] if x['in_gate'] else '-'} | {x['M2_ok'] if x['in_gate'] else '-'} |")
    print(f"\nP-M1: {msum['P_M1']}\nP-M2: {msum['P_M2']}\nreproducing WP2: {msum['n_reproducing']}/{msum['n_records']}; "
          f"resolved to 1e-6: {msum['n_resolved_1e6']}\n")
    print('### Grid convergence (A-R, P1 network, levels 0-2)\n')
    print(f"P-G1: {gsum['P_G1']}\nP-G2: {gsum['P_G2']}\nP-G3: {gsum['P_G3']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
