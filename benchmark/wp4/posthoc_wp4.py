#!/usr/bin/env python3
"""WP4 post-hoc diagnostics. None of these was pre-registered; the WP4 report labels them as post hoc.

  D1  Extended precision, what limits X2 and its G_SD counterpart at the 1e-10 level. With G(v)/G_true - 1 = a x^2 + b x^4
      (x = v/V_T) and a reference G0 = G_true (1 + delta), the Richardson residual is R = -b x1^4/100 - delta (x1 = 1 mV/V_T).
      A sinh-shaped response (diffusion-limited transport, I ~ 1 - exp(-V_D/V_T)) has a = 1/6 and b = 1/120, so
      R_sinh = -x1^4/12000. Reported: a from the +-0.1 mV value (b x^4 is 1e-12 there), R for G_DD and G_SD, R - R_sinh
      (= -delta if b = 1/120; b is not measured independently), and the WP2 forms indicator of the same mesh.
  D2  WP2 level-1 reference network values recomputed on extended-precision equilibria (P2X and AC-extended records): the
      shift relative to the WP2 P1 value, beside the WP2 rounding indicator of the same mesh (relative difference of the two
      forms, observables.same_mesh_references.*_rel_forms).
  D3  Extended precision: e_S against one unit in the last place of the largest returned current (the currents reach Python
      as float64).
  D4  DEVSIM built-in AC analysis in double precision against the WP2 P3 G_M (the DEVSIM Jacobian solved with the
      devsim.umfpack shim, phase3/c4/t1v2.umf_solve).
  D5  G_SD from the DEVSIM source contact current (Richardson) against the network G_SD, double (WP2) and extended (WP4).
  D6  Cost: wall time and peak RSS from the ledgers.
  D7  P2X_AR-4: the killed first attempt (amendment W4-A2) against the rerun, line by line.
Writes benchmark/wp4/results/posthoc_wp4.json and prints markdown tables.

Usage (repository root):  python benchmark/wp4/posthoc_wp4.py > benchmark/wp4/results/posthoc_tables_wp4.md
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PROD = os.path.join(HERE, 'results', 'production')
WP2 = os.path.join(ROOT, 'benchmark', 'wp2', 'results')
CASES = os.path.join(ROOT, 'benchmark', 'cases')
AR = ('AR-1', 'AR-2', 'AR-3', 'AR-4', 'AR-5', 'AR-6')
VT = 0.025851999786435316          # thermal voltage of the deck (300 K), as in phase3/c11/fd2d_indep.py
X1, X01 = 1e-3 / VT, 1e-4 / VT


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def ledger(path):
    out = {}
    for line in open(path):
        r = json.loads(line)
        if r['exit'] == 0:
            out[r['job']] = r
    return out


def case_level1(sid):
    for p in glob.glob(os.path.join(CASES, 'AR__*.json')):
        c = json.load(open(p))
        if c['device']['state_id'] == sid:
            smr = c['observables']['same_mesh_references']
            mid = [m for m, v in smr.items() if v['level'] == 1][0]
            return c['id'], mid, smr[mid]
    raise SystemExit(f'no A-R case for {sid}')


def rich(g1, g01):
    return (100.0 * g01 - g1) / 99.0


def gsd(d):
    I_S = {st['V_D']: st['I_source'] for st in d['states']}
    g1, g01 = -(I_S[1e-3] - I_S[-1e-3]) / 2e-3, -(I_S[1e-4] - I_S[-1e-4]) / 2e-4
    return rich(g1, g01)


def newton_iterations(d):
    """Recorded Newton iterations: equilibrium, each bias step and each return to zero bias."""
    return d['equilibrium']['iterations'] + sum(h['iterations'] for st in d['states'] for h in st['continuation']) + \
        sum(st['return_to_zero']['iterations'] for st in d['states'])


R_SINH = -X1 ** 4 / 12000.0


def d1(x):
    s = x['symmetric']
    net = x['network_same_run']
    I_S = {st['V_D']: st['I_source'] for st in x['states']}
    out = {}
    for q, g1, g01 in (('G_DD', s['1mV']['I_drain']['G'], s['0.1mV']['I_drain']['G']),
                       ('G_SD', -(I_S[1e-3] - I_S[-1e-3]) / 2e-3, -(I_S[1e-4] - I_S[-1e-4]) / 2e-4)):
        G0 = net[q]
        R = rich(g1, g01) / G0 - 1.0
        out[q] = dict(a=(g01 / G0 - 1.0) / X01 ** 2, richardson_rel=R, minus_sinh_truncation=R - R_SINH)
    out['circuit_richardson_rel'] = rich(s['1mV']['I_circuit']['G'], s['0.1mV']['I_circuit']['G']) / net['G_DD'] - 1.0
    out['weighted_richardson_rel'] = rich(s['1mV']['I_weighted']['G'], s['0.1mV']['I_weighted']['G']) / net['G_DD'] - 1.0
    return out


def d3(x):
    eS = max(abs(st['I_source'] + st['I_circuit'] + (st['I_body'] or 0.0)) for st in x['states'])
    ulp = max(float(np.spacing(max(abs(st['I_source']), abs(st['I_circuit']), abs(st['I_body'] or 0.0)))) for st in x['states'])
    eD = max(abs(st['I_drain'] - st['I_circuit']) for st in x['states'])
    return dict(e_D=eD, e_S=eS, ulp_max_current=ulp, e_S_over_ulp=eS / ulp)


def main():
    L4 = ledger(os.path.join(HERE, 'results', 'ledger.jsonl'))
    L2 = ledger(os.path.join(WP2, 'ledger.jsonl'))
    rows = []
    for sid in AR:
        cid, mid, ref = case_level1(sid)
        r = dict(state=sid, case=cid, mesh_id=mid)
        x = load(os.path.join(PROD, f'P2X_{sid}_B1_lev1.json'))
        dbl = load(os.path.join(WP2, 'production', f'P2_{sid}_B1_lev1.json'))
        p1 = load(os.path.join(WP2, 'production', f'P1_{sid}_B1_lev1.json'))
        acx = load(os.path.join(PROD, f'AC_{sid}_B1_lev1_extended.json'))
        acd = load(os.path.join(PROD, f'AC_{sid}_B1_lev1_double.json'))
        # D2: reference shift on extended equilibria
        r['D2'] = dict(G_DD_ref=ref['G_DD_S_per_cm'], G_SD_ref=ref['G_SD_S_per_cm'], rel_forms_G_DD=ref['G_DD_rel_forms'],
                       rel_forms_G_SD=ref['G_SD_rel_forms'],
                       P2X_G_DD_shift=(x['network_same_run']['G_DD'] / ref['G_DD_S_per_cm'] - 1.0) if x else None,
                       P2X_G_SD_shift=(x['network_same_run']['G_SD'] / ref['G_SD_S_per_cm'] - 1.0) if x else None,
                       ACx_G_DD_shift=(acx['network_same_run']['G_DD'] / ref['G_DD_S_per_cm'] - 1.0) if acx else None,
                       double_runs_G_DD_shift=[(d['network_same_run']['G_DD'] / ref['G_DD_S_per_cm'] - 1.0) for d in (dbl, acd) if d])
        if x:
            r['D1'] = d1(x)
            r['D3'] = d3(x)
            r['D5'] = dict(extended=gsd(x) / x['network_same_run']['G_SD'] - 1.0, double=gsd(dbl) / dbl['network_same_run']['G_SD'] - 1.0)
        if acd and p1:
            gm, gac = p1['smallsignal']['G_M'], acd['ac']['G_AC']
            r['D4'] = dict(G_AC_double=gac, G_M_WP2=gm, bitwise_equal=gac == gm, rel_diff=gac / gm - 1.0,
                           rel_G_E_WP2=p1['smallsignal']['rel_G_E_vs_G_DD'], kappa_G=p1['smallsignal']['kappa_G'],
                           E_round_estimate=p1['smallsignal']['E_round'])
        cost = {}
        for tag, job, L in (('P2_double', f'P2_{sid}_B1_lev1', L2), ('P2_extended', f'P2X_{sid}_B1_lev1', L4),
                            ('AC_double', f'AC_{sid}_B1_lev1_double', L4), ('AC_extended', f'AC_{sid}_B1_lev1_extended', L4)):
            if job in L:
                cost[tag] = dict(wall_s=L[job]['wall_s'], maxrss_MB=L[job]['maxrss_MB'])
        for tag, d in (('P2_double', dbl), ('P2_extended', x)):
            if d and tag in cost:
                cost[tag]['newton_iterations'] = newton_iterations(d)
        if acx:
            cost['AC_extended']['t_ac_solve_s'] = acx['ac']['t_ac_s']
        if acd:
            cost['AC_double']['t_ac_solve_s'] = acd['ac']['t_ac_s']
        r['D6'] = cost
        rows.append(r)
    # D7
    a1 = open(os.path.join(HERE, 'logs', 'production', 'P2X_AR-4_B1_lev1.attempt1.log')).read().splitlines()
    a2 = open(os.path.join(HERE, 'logs', 'production', 'P2X_AR-4_B1_lev1.log')).read().splitlines()
    D7 = dict(lines_attempt1=len(a1), lines_rerun=len(a2), identical_prefix=sum(1 for u, v in zip(a1, a2) if u == v),
              attempt1_is_prefix_of_rerun=a2[:len(a1)] == a1, rerun_extra_lines=a2[len(a1):])
    with open(os.path.join(HERE, 'results', 'posthoc_wp4.json'), 'w', encoding='utf-8') as f:
        json.dump(dict(note='post hoc, not pre-registered; see the docstring of benchmark/wp4/posthoc_wp4.py', V_T=VT, rows=rows, D7=D7), f, indent=1)

    w = print
    e = lambda v, p=1: 'n/a' if v is None else f'{v:+.{p}e}'  # noqa: E731
    w('### D1 (post hoc): what limits the extended-precision agreement (X2) at the 1e-10 level\n')
    w(f'R = Richardson value of the DEVSIM contact current / same-run network value − 1. a = coefficient of x² (x = V_D/V_T) from the '
      f'±0.1 mV value. A sinh-shaped response (diffusion-limited transport) has a = 1/6 = 0.16667 and leaves R_sinh = −x₁⁴/12000 = '
      f'{R_SINH:.3e}; R − R_sinh is then minus the relative error of the network value. The forms indicator is WP2\'s relative '
      'difference of the two network forms on the same mesh.\n')
    w('| State | a | R, G_DD | R, G_SD | R − R_sinh, G_DD | R − R_sinh, G_SD | WP2 forms indicator G_DD / G_SD | R, circuit / weighting-function current |')
    w('|---|---:|---:|---:|---:|---:|---|---|')
    for r in rows:
        if 'D1' in r:
            q, f2 = r['D1'], r['D2']
            w(f"| {r['state']} | {q['G_DD']['a']:.5f} | {q['G_DD']['richardson_rel']:+.2e} | {q['G_SD']['richardson_rel']:+.2e} | "
              f"{q['G_DD']['minus_sinh_truncation']:+.1e} | {q['G_SD']['minus_sinh_truncation']:+.1e} | "
              f"{f2['rel_forms_G_DD']:+.1e} / {f2['rel_forms_G_SD']:+.1e} | {q['circuit_richardson_rel']:+.2e} / {q['weighted_richardson_rel']:+.2e} |")
    w('\n### D2 (post hoc): WP2 level-1 reference values recomputed on extended-precision equilibria\n')
    w('| State | G_DD shift (P2X) | G_DD shift (AC ext.) | WP2 forms indicator G_DD | G_SD shift (P2X) | WP2 forms indicator G_SD | double-precision reruns |')
    w('|---|---:|---:|---:|---:|---:|---|')
    for r in rows:
        q = r['D2']
        w(f"| {r['state']} | {e(q['P2X_G_DD_shift'])} | {e(q['ACx_G_DD_shift'])} | {e(q['rel_forms_G_DD'])} | {e(q['P2X_G_SD_shift'])} | "
          f"{e(q['rel_forms_G_SD'])} | {', '.join(e(v) for v in q['double_runs_G_DD_shift'])} |")
    w('\n### D3 (post hoc): extended-precision current errors against float64 resolution\n')
    w('| State | e_D (A/cm) | e_S (A/cm) | ulp of largest current (A/cm) | e_S / ulp |')
    w('|---|---:|---:|---:|---:|')
    for r in rows:
        if 'D3' in r:
            q = r['D3']
            w(f"| {r['state']} | {q['e_D']:.1e} | {q['e_S']:.1e} | {q['ulp_max_current']:.1e} | {q['e_S_over_ulp']:.2f} |")
    w('\n### D4 (post hoc): DEVSIM built-in AC (double) against the WP2 small-signal `G_M` (same Jacobian, UMFPACK through `devsim.umfpack`)\n')
    w('| State | G_AC double (S/cm) | WP2 G_M (S/cm) | bitwise equal | relative difference | WP2 κ_G | WP2 E_round estimate |')
    w('|---|---:|---:|---|---:|---:|---:|')
    for r in rows:
        if 'D4' in r:
            q = r['D4']
            w(f"| {r['state']} | {q['G_AC_double']:.6e} | {q['G_M_WP2']:.6e} | {'yes' if q['bitwise_equal'] else '**no**'} | {q['rel_diff']:+.1e} | "
              f"{q['kappa_G']:.1e} | {q['E_round_estimate']:.1e} |")
    w('\n### D5 (post hoc): G_SD from the DEVSIM source contact current (Richardson) vs the same-run network G_SD\n')
    w('| State | double (WP2 P2) | extended (WP4 P2X) |')
    w('|---|---:|---:|')
    for r in rows:
        if 'D5' in r:
            w(f"| {r['state']} | {r['D5']['double']:+.2e} | {r['D5']['extended']:+.2e} |")
    w('\n### D6 (post hoc): cost (ledger wall time and peak RSS)\n')
    w('| State | P2 double | P2 extended | ratio | Newton iterations double / extended | AC double (solve) | AC extended (solve) |')
    w('|---|---:|---:|---:|---:|---:|---:|')
    for r in rows:
        c = r['D6']
        f = lambda t: f"{c[t]['wall_s']:.0f} s, {c[t]['maxrss_MB']} MB" if t in c else 'n/a'  # noqa: E731
        g = lambda t: f"{c[t]['wall_s']:.0f} s ({c[t]['t_ac_solve_s']:.1f} s), {c[t]['maxrss_MB']} MB" if t in c else 'n/a'  # noqa: E731
        ratio = f"{c['P2_extended']['wall_s'] / c['P2_double']['wall_s']:.1f}×" if 'P2_extended' in c and 'P2_double' in c else 'n/a'
        its = (f"{c['P2_double']['newton_iterations']} / {c['P2_extended']['newton_iterations']}"
               if 'P2_extended' in c and 'P2_double' in c else 'n/a')
        w(f"| {r['state']} | {f('P2_double')} | {f('P2_extended')} | {ratio} | {its} | {g('AC_double')} | {g('AC_extended')} |")
    w(f"\n### D7 (post hoc): P2X_AR-4, killed first attempt vs rerun\n\nThe first attempt's log has {D7['lines_attempt1']} lines; "
      f"{D7['identical_prefix']} of them are identical to the rerun's first {D7['lines_attempt1']} lines "
      f"(attempt 1 is {'a' if D7['attempt1_is_prefix_of_rerun'] else '**not** a'} prefix of the rerun). The rerun adds {len(D7['rerun_extra_lines'])} "
      'lines: the summary, the output path and the exit code.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
