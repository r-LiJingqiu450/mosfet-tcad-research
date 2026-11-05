#!/usr/bin/env python3
"""Sensitivity of the WP6 results to the audit amendment W6-A2 and to contamination (pre-registration section 8).

1. Original audit: runs audit.py as frozen (commit 710e55a, unchanged since 01195f2) on every run of the stage. It writes
   audit_v1.json where missing and checks that it reproduces a stored audit_v1.json where present (W6-001 to W6-011 were
   audited with it before the amendment). The frozen code is loaded from git, so nothing in the repository is changed;
   it writes into a scratch copy of the run record and never touches the run's audit.json.
2. Evaluation of the stage with the original audit instead of the amended one (same classification code).
3. Registered worst case: contaminated runs counted as SE in their arm, under the amended and under the original audit.

Writes results/sensitivity_wp6_stage<k>.json and prints a Markdown table.

Usage (repository root):  python benchmark/wp6/harness/sensitivity_wp6.py --stage 1
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
WP6 = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(WP6))
RES = os.path.join(WP6, 'results')
RUNS = os.path.join(RES, 'runs')
FROZEN = '710e55a'
sys.path.insert(0, HERE)
import evaluate_wp6 as EV  # noqa: E402


def frozen_audit_module():
    src = subprocess.check_output(['git', 'show', f'{FROZEN}:benchmark/wp6/harness/audit.py'], cwd=ROOT, text=True)
    mod = types.ModuleType('audit_frozen')
    mod.__file__ = os.path.join(HERE, 'audit.py')      # so that ROOT and the guard import resolve as at the freeze
    exec(compile(src, f'audit.py@{FROZEN}', 'exec'), mod.__dict__)
    return mod


def audit_v1(mod, rd):
    """Run the frozen audit on a scratch copy of the run record; store or verify audit_v1.json."""
    with tempfile.TemporaryDirectory() as tmp:
        for f in ('run.json', 'checker_log.jsonl'):
            if os.path.exists(os.path.join(rd, f)):
                shutil.copy2(os.path.join(rd, f), tmp)
        a = mod.audit(tmp)
    p = os.path.join(rd, 'audit_v1.json')
    if os.path.exists(p):
        stored = json.load(open(p))
        same = json.loads(json.dumps(a)) == stored
        return stored, ('reproduced' if same else 'DIFFERS from stored')
    with open(p, 'w') as f:
        json.dump(a, f, indent=1)
    return a, 'written'


def rows_with(audits, stage):
    tasks = {t['task']: t for t in json.load(open(os.path.join(WP6, 'tasks', 'tasks.json')))['tasks']}
    rows = []
    for rid, au in sorted(audits.items()):
        run = json.load(open(os.path.join(RUNS, rid, 'run.json')))
        if run['stage'] != stage:
            continue
        t = tasks[run['task']]
        r = EV.classify(run, au, t['reference_G_DD_S_per_cm'])
        r['task_class'] = 'low' if t['reference_G_DD_S_per_cm'] * 1e-3 <= EV.LOW_CURRENT_A else 'moderate'
        rows.append(r)
    return rows


def worst_case(rows):
    out = copy.deepcopy(rows)
    for r in out:
        if r['exclusion'] and r['exclusion'].startswith('contaminated'):
            r['valid'], r['outcome'], r['worst_case_from'] = True, 'SE', r['exclusion']
    return out


def brief(rows):
    s = EV.summarize(rows)
    keep = {}
    for arm in ('A', 'B'):
        x = s[arm]
        if x['n_runs']:
            keep[arm] = dict(n_runs=x['n_runs'], n_valid=x['n_valid'], n_contaminated=x['n_contaminated'], counts=x['counts'],
                             SE_rate=x['SE_rate'], SE_CI95=x['SE_CI95'])
    keep['stop_rule_stage1'] = s['stop_rule_stage1']['applies']
    keep['halt_rule_stage1'] = s['halt_rule_stage1']
    keep['halt_rule_stage2'] = s['halt_rule_stage2']
    keep['H4'] = dict(p_value=s['H4']['p_value'], supported=s['H4']['supported'])
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', type=int, required=True)
    a = ap.parse_args()
    mod = frozen_audit_module()
    amended, original, checks = {}, {}, {}
    for rid in sorted(os.listdir(RUNS)):
        rd = os.path.join(RUNS, rid)
        if not os.path.exists(os.path.join(rd, 'run.json')) or not rid.startswith('W6-'):
            continue
        if json.load(open(os.path.join(rd, 'run.json')))['stage'] != a.stage:
            continue
        amended[rid] = json.load(open(os.path.join(rd, 'audit.json')))
        original[rid], checks[rid] = audit_v1(mod, rd)
    variants = {
        'primary (amended audit W6-A2)': rows_with(amended, a.stage),
        'original audit (frozen 710e55a)': rows_with(original, a.stage),
    }
    variants['worst case, amended audit: contaminated = SE'] = worst_case(variants['primary (amended audit W6-A2)'])
    variants['worst case, original audit: contaminated = SE'] = worst_case(variants['original audit (frozen 710e55a)'])
    changed = [dict(run_id=rid, amended=dict(contaminated=amended[rid]['contaminated'], F1=amended[rid]['flags']['F1_forbidden_access']),
                    original=dict(contaminated=original[rid]['contaminated'], F1=original[rid]['flags']['F1_forbidden_access']))
               for rid in amended if amended[rid]['contaminated'] != original[rid]['contaminated']]
    out = dict(stage=a.stage, frozen_audit_commit=FROZEN, audit_v1_status=checks, runs_whose_validity_differs=changed,
               variants={k: brief(v) for k, v in variants.items()},
               worst_case_rows={k: [r for r in v if r.get('worst_case_from')] for k, v in variants.items() if k.startswith('worst')})
    with open(os.path.join(RES, f'sensitivity_wp6_stage{a.stage}.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(f'### Sensitivity, stage {a.stage}\n')
    print('| Variant | arm | valid | SE | CS | FE | FC | NR | SE rate (95% CI) | stop rule applies |')
    print('|---|---|---:|---:|---:|---:|---:|---:|---|---|')
    for k, v in out['variants'].items():
        for arm in ('A', 'B'):
            if arm not in v:
                continue
            x = v[arm]
            c, ci = x['counts'], x['SE_CI95']
            rate = 'n/a' if x['SE_rate'] is None else f"{x['SE_rate']:.3f} ({ci[0]:.3f}–{ci[1]:.3f})"
            print(f"| {k} | {arm} | {x['n_valid']} | {c['SE']} | {c['CS']} | {c['FE']} | {c['FC']} | {c['NR']} | {rate} | "
                  f"{v['stop_rule_stage1'] if a.stage == 1 else 'n/a'} |")
    print(f"\nRuns whose validity differs between the audits: {[c['run_id'] for c in changed]}")
    print(f"audit_v1.json: {sum(s == 'reproduced' for s in checks.values())} reproduced, "
          f"{sum(s == 'written' for s in checks.values())} written, {[k for k, s in checks.items() if s.startswith('DIFF')]} differ")
    return 0


if __name__ == '__main__':
    sys.exit(main())
