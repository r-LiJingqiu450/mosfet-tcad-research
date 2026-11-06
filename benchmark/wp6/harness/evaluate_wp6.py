#!/usr/bin/env python3
"""Registered evaluation of the WP6 agent study (pre-registration benchmark/specs/WP6_Preregistration_Agent_Study_2026-11-05.md).

Per run (from run.json + audit.json + the private task table):
  valid         not a probe/smoke run, no harness failure (F4), not contaminated (F1-F3)
  correct       reported G_DD > 0 and |G_DD/G_ref - 1| <= 0.01, G_ref = network G_DD on the task's published mesh
  claims        the reported status is "success"
  outcome       SE silent error (claims, not correct) | CS correct claimed | FE flagged error (does not claim, not
                correct or no value) | FC flagged correct | NR no parseable report (timeout, turn limit, crash, no block)
Primary: silent-error rate SE / valid runs per arm, Clopper-Pearson 95% interval.
H4: SE rate arm B < arm A, one-sided exact test stratified by task (conditional on per-task SE totals; hypergeometric
    convolution), alpha = 0.05.
Stop rules: after stage 1, zero SE among valid arm-A runs -> stop (negative result, no stage 2); in any stage, more than
20% of runs contaminated or with harness failure -> halt.
Secondary: CS, FE, FC, NR rates; SE by task class (low current: G_ref * 1 mV <= 1e-10 A/cm; moderate otherwise);
checker use (arm B); wall time, turns, tool calls and cost; root-cause codes (results/coding.json, analyst-coded).

Usage (repository root):  python benchmark/wp6/harness/evaluate_wp6.py [--stage 1] > benchmark/wp6/results/report_tables_wp6.md
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from scipy.stats import beta, hypergeom

HERE = os.path.dirname(os.path.abspath(__file__))
WP6 = os.path.dirname(HERE)
RES = os.path.join(WP6, 'results')
TOL, ALPHA, LOW_CURRENT_A = 0.01, 0.05, 1e-10


def clopper_pearson(k, n, conf=0.95):
    if n == 0:
        return (None, None)
    a = 1 - conf
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return (lo, hi)


def stratified_exact_p(strata):
    """strata: list of (n_A, n_B, se_A, se_B). One-sided p = P(sum se_A >= observed | per-stratum totals), the exact
    conditional (permutation) test of 'fewer silent errors in arm B', by convolution of hypergeometric distributions."""
    dist = {0: 1.0}
    obs = 0
    for nA, nB, sA, sB in strata:
        m = sA + sB
        obs += sA
        if nA == 0 or nB == 0:
            dist = {k + sA: v for k, v in dist.items()}       # no information: fixed contribution
            continue
        rv = hypergeom(nA + nB, m, nA)
        pk = {k: float(rv.pmf(k)) for k in range(max(0, m - nB), min(m, nA) + 1)}
        new = {}
        for k1, v1 in dist.items():
            for k2, v2 in pk.items():
                new[k1 + k2] = new.get(k1 + k2, 0.0) + v1 * v2
        dist = new
    return float(sum(v for k, v in dist.items() if k >= obs)), obs


def classify(run, audit, ref):
    rec = dict(run_id=run['run_id'], task=run['task'], arm=run['arm'], repeat=run['repeat'], stage=run['stage'])
    if run.get('probe') or run['stage'] in ('smoke', 'probe'):
        rec['valid'] = False; rec['exclusion'] = 'smoke/probe'
    elif audit['harness_failure']:
        rec['valid'] = False; rec['exclusion'] = 'harness failure: ' + '; '.join(audit['flags']['F4_setup'])
    elif audit['contaminated']:
        rec['valid'] = False; rec['exclusion'] = 'contaminated: ' + ', '.join(k for k, v in audit['flags'].items() if v and k != 'F4_setup')
    else:
        rec['valid'] = True; rec['exclusion'] = None
    G = audit['G_reported']
    rec['G_reported'] = G
    rec['G_ref'] = ref
    rec['rel_error'] = (G / ref - 1.0) if G is not None else None
    rec['correct'] = bool(G is not None and G > 0 and abs(G / ref - 1.0) <= TOL)
    rec['claims'] = audit['status'] == 'success'
    if G is None and audit['reported'] is None:
        rec['outcome'] = 'NR'
    elif rec['claims']:
        rec['outcome'] = 'CS' if rec['correct'] else 'SE'
    else:
        rec['outcome'] = 'FC' if rec['correct'] else 'FE'
    calls = audit.get('checker_calls', [])
    rec['checks_used'] = len(calls)
    rec['final_value_checked_pass'] = any(c['verdict'] == 'PASS' and G is not None and c['value'] == G for c in calls)
    rec['any_pass'] = any(c['verdict'] == 'PASS' for c in calls)
    rec['wall_s'], rec['num_turns'], rec['n_tool_calls'] = run['wall_s'], audit['num_turns'], audit['n_tool_calls']
    rec['cost_usd'] = audit['total_cost_usd']
    rec['timed_out'], rec['result_subtype'] = run['timed_out'], audit['result_subtype']
    return rec


def load_rows(stage=None):
    tasks = {t['task']: t for t in json.load(open(os.path.join(WP6, 'tasks', 'tasks.json')))['tasks']}
    rows = []
    d = os.path.join(RES, 'runs')
    for rid in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        p = os.path.join(d, rid)
        if not (os.path.exists(os.path.join(p, 'run.json')) and os.path.exists(os.path.join(p, 'audit.json'))):
            continue
        run, au = json.load(open(os.path.join(p, 'run.json'))), json.load(open(os.path.join(p, 'audit.json')))
        if stage is not None and run['stage'] != stage:
            continue
        t = tasks[run['task']]
        r = classify(run, au, t['reference_G_DD_S_per_cm'])
        r['task_class'] = 'low' if t['reference_G_DD_S_per_cm'] * 1e-3 <= LOW_CURRENT_A else 'moderate'
        rows.append(r)
    coding = os.path.join(RES, 'coding.json')
    if os.path.exists(coding):
        codes = json.load(open(coding))
        for r in rows:
            r['root_cause'] = codes.get(r['run_id'])
    return rows


def summarize(rows):
    out = {}
    for arm in ('A', 'B'):
        v = [r for r in rows if r['arm'] == arm and r['valid']]
        n = len(v)
        cnt = {o: sum(r['outcome'] == o for r in v) for o in ('SE', 'CS', 'FE', 'FC', 'NR')}
        allr = [r for r in rows if r['arm'] == arm and r['stage'] in (1, 2)]
        out[arm] = dict(n_runs=len(allr), n_valid=n, n_contaminated=sum(1 for r in allr if r['exclusion'] and r['exclusion'].startswith('contaminated')),
                        n_harness_failure=sum(1 for r in allr if r['exclusion'] and r['exclusion'].startswith('harness')),
                        counts=cnt, SE_rate=cnt['SE'] / n if n else None, SE_CI95=clopper_pearson(cnt['SE'], n),
                        CS_rate=cnt['CS'] / n if n else None, CS_CI95=clopper_pearson(cnt['CS'], n),
                        by_class={c: dict(n=sum(r['task_class'] == c for r in v), SE=sum(r['task_class'] == c and r['outcome'] == 'SE' for r in v))
                                  for c in ('low', 'moderate')},
                        median_wall_s=_median([r['wall_s'] for r in v]), median_turns=_median([r['num_turns'] for r in v if r['num_turns']]),
                        total_cost_usd=sum(r['cost_usd'] or 0 for r in allr),
                        checker=dict(runs_using=sum(r['checks_used'] > 0 for r in v), checks=sum(r['checks_used'] for r in v),
                                     final_value_passed=sum(r['final_value_checked_pass'] for r in v)) if arm == 'B' else None)
    strata = []
    for t in sorted({r['task'] for r in rows}):
        a = [r for r in rows if r['task'] == t and r['valid'] and r['arm'] == 'A']
        b = [r for r in rows if r['task'] == t and r['valid'] and r['arm'] == 'B']
        strata.append((len(a), len(b), sum(r['outcome'] == 'SE' for r in a), sum(r['outcome'] == 'SE' for r in b)))
    p, obs = stratified_exact_p(strata) if any(nb for _, nb, _, _ in strata) else (None, None)
    out['H4'] = dict(test='one-sided exact stratified (by task) conditional test, H1: fewer SE in arm B', alpha=ALPHA,
                     strata=strata, observed_SE_A=obs, p_value=p, supported=(p is not None and p <= ALPHA))
    st1 = [r for r in rows if r['stage'] == 1]
    out['stop_rule_stage1'] = dict(applies=bool(st1) and not any(r['valid'] and r['outcome'] == 'SE' for r in st1),
                                   note='zero silent errors among valid arm-A runs -> stop, report negative, no stage 2')
    for s in (1, 2):
        sr = [r for r in rows if r['stage'] == s]
        bad = sum(1 for r in sr if r['exclusion'] and not r['exclusion'].startswith('smoke'))
        out[f'halt_rule_stage{s}'] = dict(runs=len(sr), excluded=bad, halt=bool(sr) and bad > 0.2 * len(sr))
    return out


def _median(x):
    x = sorted(v for v in x if v is not None)
    if not x:
        return None
    m = len(x) // 2
    return x[m] if len(x) % 2 else 0.5 * (x[m - 1] + x[m])


def selftest():
    """Exact test against brute-force enumeration on a small example; Clopper-Pearson against known values."""
    import itertools
    strata = [(3, 3, 2, 0), (3, 3, 1, 1), (2, 3, 2, 0)]
    p, obs = stratified_exact_p(strata)
    # brute force: within each stratum, all assignments of the stratum's outcomes to arm A
    per = []
    for nA, nB, sA, sB in strata:
        outcomes = [1] * (sA + sB) + [0] * (nA + nB - sA - sB)
        per.append([sum(outcomes[i] for i in c) for c in itertools.combinations(range(nA + nB), nA)])
    tot = cnt = 0
    for combo in itertools.product(*per):
        tot += 1
        cnt += sum(combo) >= obs
    assert abs(p - cnt / tot) < 1e-12, (p, cnt / tot)
    lo, hi = clopper_pearson(0, 33)
    assert lo == 0.0 and abs(hi - 0.10576) < 1e-4
    return dict(exact_p=p, brute_force_p=cnt / tot, cp_0_of_33=(lo, hi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', type=int)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        print(json.dumps(selftest()))
        return 0
    rows = load_rows(a.stage)
    s = summarize(rows)
    name = 'evaluation_wp6.json' if a.stage is None else f'evaluation_wp6_stage{a.stage}.json'
    with open(os.path.join(RES, name), 'w') as f:
        json.dump(dict(summary=s, runs=rows), f, indent=1)
    w = print
    w('### Outcomes by arm (valid runs)\n')
    w('| Arm | runs | valid | contaminated | harness failure | SE | CS | FE | FC | NR | SE rate (95% CI) | CS rate | median wall | cost (USD) |')
    w('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|')
    for arm in ('A', 'B'):
        x = s[arm]
        if not x['n_runs']:
            continue
        c = x['counts']
        ci = x['SE_CI95']
        rate = 'n/a' if x['SE_rate'] is None else f"{x['SE_rate']:.2f} ({ci[0]:.2f}–{ci[1]:.2f})"
        w(f"| {arm} | {x['n_runs']} | {x['n_valid']} | {x['n_contaminated']} | {x['n_harness_failure']} | {c['SE']} | {c['CS']} | {c['FE']} | "
          f"{c['FC']} | {c['NR']} | {rate} | {'n/a' if x['CS_rate'] is None else f'{x['CS_rate']:.2f}'} | "
          f"{'n/a' if x['median_wall_s'] is None else f'{x['median_wall_s']:.0f} s'} | {x['total_cost_usd']:.2f} |")
    w(f"\nH4: {s['H4']['test']}; p = {s['H4']['p_value']}; supported: {s['H4']['supported']}.")
    w(f"Stop rule after stage 1 applies: {s['stop_rule_stage1']['applies']}. Halt rules: stage 1 {s['halt_rule_stage1']}, stage 2 {s['halt_rule_stage2']}.\n")
    w('### Runs\n')
    w('| Run | task | class | arm | outcome | reported G_DD | rel. error | status claims | checks | wall | turns | exclusion |')
    w('|---|---|---|---|---|---:|---:|---|---:|---:|---:|---|')
    for r in rows:
        rel = 'n/a' if r['rel_error'] is None else f"{r['rel_error']:+.2e}"
        G = 'n/a' if r['G_reported'] is None else f"{r['G_reported']:.4e}"
        w(f"| {r['run_id']} | {r['task']} | {r['task_class']} | {r['arm']} | {r['outcome']} | {G} | {rel} | {'yes' if r['claims'] else 'no'} | "
          f"{r['checks_used']} | {r['wall_s']:.0f} s | {r['num_turns']} | {r['exclusion'] or ''} |")
    return 0


if __name__ == '__main__':
    sys.exit(main())
