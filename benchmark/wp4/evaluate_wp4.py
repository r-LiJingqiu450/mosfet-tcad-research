#!/usr/bin/env python3
"""Evaluate WP4 against the frozen pre-registration (benchmark/specs/WP4_Preregistration_External_Paths_2026-11-04.md).

  WP4-A  H3: the DEVSIM default drain contact current (WP2 P2 Richardson values) scored with tools/score.py on the A-R, HO
         and B-R level-1 cases (formal record; outcome known from the WP2 examples).
  X1     extended precision: max over biases of e_D = |I_D - I_C| and e_S = |I_S + I_C + I_B| <= 1e-20 A/cm, per run.
  X2     extended precision: Richardson G_DD from the DEVSIM default drain current within 1e-6 of the same-run network.
  X3     extended precision: circuit-current Richardson within 1e-6 of the same-run network (E-CIRC).
  Y1/Y2  DEVSIM built-in AC analysis, double / extended: |G_AC/G_DD - 1| > 1e-2 for AR-1, AR-2, AR-3.
Each is reported per run as n/6 (X) or n/3 (Y); none is a kill. Writes benchmark/wp4/results/evaluation_wp4.json and
prints markdown tables (benchmark/wp4/results/report_tables_wp4.md when redirected).
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PROD = os.path.join(HERE, 'results', 'production')
WP2 = os.path.join(ROOT, 'benchmark', 'wp2', 'results', 'production')
sys.path.insert(0, os.path.join(ROOT, 'benchmark', 'tools'))
import score as S  # noqa: E402

AR = ('AR-1', 'AR-2', 'AR-3', 'AR-4', 'AR-5', 'AR-6')
PSI = dict(zip(AR, (12, 16, 21, 12, 16, 21)))
TAU_NUM, LIM_ABS, LIM_AC = 1e-6, 1e-20, 1e-2


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def errors(d):
    st = d['states']
    return (max(abs(s['I_drain'] - s['I_circuit']) for s in st),
            max(abs(s['I_source'] + s['I_circuit'] + (s['I_body'] or 0.0)) for s in st))


def rich(d, k):
    s = d['symmetric']
    return (100.0 * s['0.1mV'][k]['G'] - s['1mV'][k]['G']) / 99.0


def b1_rows():
    rows = []
    for sid in AR:
        x = load(os.path.join(PROD, f'P2X_{sid}_B1_lev1.json'))
        dbl = load(os.path.join(WP2, f'P2_{sid}_B1_lev1.json'))
        if x is None:
            rows.append(dict(state=sid, missing=True))
            continue
        G = x['network_same_run']['G_DD']
        eD, eS = errors(x)
        eD0, eS0 = errors(dbl)
        r = dict(state=sid, psi_hat=PSI[sid], G_DD=G, G_DD_double_run=dbl['network_same_run']['G_DD'],
                 network_ext_vs_double=G / dbl['network_same_run']['G_DD'] - 1.0,
                 e_D=eD, e_S=eS, e_D_double=eD0, e_S_double=eS0, reduction_D=eD0 / eD if eD > 0 else float('inf'),
                 default_1mV=x['symmetric']['1mV']['I_drain']['G'] / G - 1.0, default_0p1mV=x['symmetric']['0.1mV']['I_drain']['G'] / G - 1.0,
                 default_richardson=rich(x, 'I_drain') / G - 1.0, default_richardson_double=rich(dbl, 'I_drain') / dbl['network_same_run']['G_DD'] - 1.0,
                 circuit_richardson=rich(x, 'I_circuit') / G - 1.0, weighted_richardson=rich(x, 'I_weighted') / G - 1.0,
                 newton_final_update_max=max(max(h['last'][-1]['dpsi'], h['last'][-1]['dphin']) for s in x['states'] for h in s['continuation']),
                 wall_s=x['wall_s'], wall_s_double=dbl['wall_s'],
                 flags=x['precision']['parameters_after_run'])
        r['X1'] = eD <= LIM_ABS and eS <= LIM_ABS
        r['X2'] = abs(r['default_richardson']) <= TAU_NUM
        r['X3'] = abs(r['circuit_richardson']) <= TAU_NUM
        rows.append(r)
    done = [r for r in rows if not r.get('missing')]
    summ = {k: f"{sum(r[k] for r in done)}/{len(done)}" for k in ('X1', 'X2', 'X3')}
    return rows, summ


def b2_rows():
    rows = []
    for sid in AR:
        r = dict(state=sid, psi_hat=PSI[sid])
        for prec in ('double', 'extended'):
            d = load(os.path.join(PROD, f'AC_{sid}_B1_lev1_{prec}.json'))
            if d is None:
                r[prec] = None
                continue
            r[prec] = dict(G_AC=d['ac']['G_AC'], rel=d['ac']['rel_vs_network'], I_imag=d['ac']['I_imag'], wall_s=d['wall_s'],
                           G_DD=d['network_same_run']['G_DD'])
        p3 = load(os.path.join(WP2, f'P1_{sid}_B1_lev1.json'))
        r['P3_G_M_rel_WP2'] = p3['smallsignal']['rel_G_M_vs_G_DD'] if p3 else None
        r['P3_G_E_rel_WP2'] = p3['smallsignal']['rel_G_E_vs_G_DD'] if p3 else None
        r['predicted'] = sid in ('AR-1', 'AR-2', 'AR-3')
        for k, prec in (('Y1', 'double'), ('Y2', 'extended')):
            r[k] = (abs(r[prec]['rel']) > LIM_AC) if (r['predicted'] and r[prec]) else None
        rows.append(r)
    pred = [r for r in rows if r['predicted']]
    summ = {k: f"{sum(bool(r[k]) for r in pred if r[k] is not None)}/{sum(r[k] is not None for r in pred)}" for k in ('Y1', 'Y2')}
    return rows, summ


def h3():
    cases = S.load_cases()
    res = {}
    for k, c in sorted(cases.items()):
        if c['family'] == 'zero_bias_3terminal':
            lev = c['mesh']['level']
            e = [x for x in c['known_bad_outputs'] + c['cross_paths'] if x['name'] == f'devsim_drain_contact_current_richardson_lev{lev}']
            if e:
                res[k] = dict(G_DD_S_per_cm=e[0]['G_DD_S_per_cm'], same_mesh=True)
        elif c['tier'] == 'B-R' and c['mesh']['level'] == 1:
            e = [x for x in c['known_bad_outputs'] + c['cross_paths'] if x['name'] == 'devsim_drain_contact_current_richardson_lev1']
            if e:
                res[k] = dict(G0_S_per_cm=e[0]['G0_S_per_cm'], same_mesh=True)
    sub = dict(submission=dict(name='DEVSIM 2.11.2 default drain contact current (get_contact_current), double precision',
                               current_evaluation='symmetric +-0.1/+-1 mV differences, Richardson'), results=res)
    rep = S.score(sub, cases)
    by_tier = {}
    for r in rep['results']:
        t = by_tier.setdefault(r['tier'], {})
        t[r['status']] = t.get(r['status'], 0) + 1
    trivial = all(r['status'] == 'pass' for r in rep['results'])
    return dict(submission=sub, report=rep, by_tier=by_tier, kill_every_case_passes=trivial)


def e(x, p=2):
    return 'n/a' if x is None else f'{x:+.{p}e}'


def main():
    B1, s1 = b1_rows()
    B2, s2 = b2_rows()
    H3 = h3()
    out = dict(B1=B1, B1_summary=s1, B2=B2, B2_summary=s2, H3=dict(by_tier=H3['by_tier'], kill=H3['kill_every_case_passes'],
               results=H3['report']['results']))
    with open(os.path.join(HERE, 'results', 'evaluation_wp4.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(HERE, 'results', 'h3_devsim_default_submission.json'), 'w', encoding='utf-8') as f:
        json.dump(H3['submission'], f, indent=1, sort_keys=True)
    w = print
    w('### WP4-A: H3, DEVSIM default drain contact current scored with `tools/score.py`\n')
    w('| Tier | result |')
    w('|---|---|')
    for t, c in sorted(H3['by_tier'].items()):
        w(f"| {t} | {', '.join(f'{k} {v}' for k, v in sorted(c.items()))} |")
    w(f"\nH3 kill (every external-path case passes trivially): {H3['kill_every_case_passes']}.\n")
    w('### WP4-B1: DEVSIM extended precision, post-solve contact current (level 1, body B1; paired with the WP2 double-precision run)\n')
    w('| State | ψ̂ | e_D double (A/cm) | e_D extended | e_S extended | reduction | DEVSIM drain G, Richardson: double → extended | circuit Richardson (ext) | X1 | X2 | X3 | wall double → extended |')
    w('|---|---:|---:|---:|---:|---:|---|---:|---|---|---|---|')
    for r in B1:
        if r.get('missing'):
            w(f"| {r['state']} | missing |||||||||||")
            continue
        yn = lambda v: 'yes' if v else '**no**'  # noqa: E731
        w(f"| {r['state']} | {r['psi_hat']} | {r['e_D_double']:.1e} | {r['e_D']:.1e} | {r['e_S']:.1e} | {r['reduction_D']:.1e} | "
          f"{e(r['default_richardson_double'])} → {e(r['default_richardson'])} | {e(r['circuit_richardson'], 1)} | {yn(r['X1'])} | {yn(r['X2'])} | "
          f"{yn(r['X3'])} | {r['wall_s_double']:.0f} s → {r['wall_s']:.0f} s |")
    w(f"\nB1 summary: X1 {s1['X1']}, X2 {s1['X2']}, X3 {s1['X3']}.")
    done = [r for r in B1 if not r.get('missing')]
    if done:
        w(f"Network G_DD from the extended equilibrium vs the double equilibrium: |rel| ≤ {max(abs(r['network_ext_vs_double']) for r in done):.1e}; "
          f"final Newton updates ≤ {max(r['newton_final_update_max'] for r in done):.1e} V.\n")
    w('### WP4-B2: DEVSIM built-in AC analysis at 1e-6 Hz (level 1, body B1), G_AC vs the same-run network G_DD\n')
    w('| State | ψ̂ | double | extended | WP2 P3 `G_M` (UMFPACK) | predicted > 1e-2 | Y1 (double) | Y2 (extended) |')
    w('|---|---:|---:|---:|---:|---|---|---|')
    for r in B2:
        dd = r['double']['rel'] if r['double'] else None
        xx = r['extended']['rel'] if r['extended'] else None
        f = lambda v: '—' if v is None else ('confirmed' if v else '**refuted**')  # noqa: E731
        w(f"| {r['state']} | {r['psi_hat']} | {e(dd)} | {e(xx)} | {e(r['P3_G_M_rel_WP2'])} | {'yes' if r['predicted'] else 'no'} | {f(r['Y1'])} | {f(r['Y2'])} |")
    w(f"\nB2 summary: Y1 {s2['Y1']}, Y2 {s2['Y2']}.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
