#!/usr/bin/env python3
"""Exploratory (post hoc, not pre-registered) description of a WP6 stage: what the agents did. Decides no outcome.

From the run records, audits, transcripts and the analyst's method coding (results/method_coding_stage<k>.json):
  1. per task: class, G_ref, outcomes, largest |relative error| of the reported values, wall times;
  2. reported values: distribution of the relative error;
  3. paths the agents' code touched (audit.method_hints, counted from tool inputs; a mechanical count);
  4. final method, by category (analyst coding);
  5. refused tool calls, by the mechanism that refused them, read from the tool results in the transcripts:
     the workspace guard (its block message), Claude Code's own command check (<tool_use_error>), and the auto-mode
     permission classifier (permission_denied events).
Writes results/exploratory_wp6_stage<k>.json and prints Markdown tables.

Usage (repository root):  python benchmark/wp6/harness/explore_wp6.py --stage 1
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WP6 = os.path.dirname(HERE)
RES = os.path.join(WP6, 'results')
RUNS = os.path.join(RES, 'runs')
sys.path.insert(0, HERE)
import evaluate_wp6 as EV  # noqa: E402

GUARD_MARK = 'Blocked by workspace policy'
HINT_LABELS = {'get_contact_current': 'DEVSIM contact current (get_contact_current)', 'ac_analysis': "AC small-signal solve (type='ac')",
               'circuit_node_value': 'circuit current V1.I (get_circuit_node_value)', 'circuit_alter': 'drain bias steps (circuit_alter)',
               'extended_precision': 'extended precision (extended_solver/model/equation)', 'jacobian_export': 'Jacobian export (get_matrix_and_rhs)',
               'richardson_or_extrapolation': 'Richardson / extrapolation (by name)', 'electron_density_postprocessing': 'electron density read back (Electrons)'}


def refusals(run_id):
    """Refused tool calls of one run, by mechanism, from the transcript."""
    p = os.path.join(RUNS, run_id, 'transcript.jsonl.gz')
    out = collections.Counter()
    targets = collections.Counter()
    with gzip.open(p, 'rt') as f:
        for line in f:
            e = json.loads(line)
            if e.get('type') == 'system' and e.get('subtype') == 'permission_denied':
                out[f"classifier: {e.get('decision_reason')}"] += 1
            if e.get('type') != 'user':
                continue
            for c in e.get('message', {}).get('content', []) or []:
                if not (isinstance(c, dict) and c.get('type') == 'tool_result' and c.get('is_error')):
                    continue
                t = c.get('content')
                t = t if isinstance(t, str) else ' '.join(x.get('text', '') for x in t if isinstance(x, dict))
                if GUARD_MARK in t:
                    out['workspace guard'] += 1
                    targets['own background-task output' if '/tmp/claude-' in t or 'temporary directory' in t
                            else ('code text matched the home-directory rule' if 'home directory' in t else 'other')] += 1
                elif t.startswith('<tool_use_error>'):
                    out["Claude Code's command check"] += 1
    return dict(out), dict(targets)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', type=int, required=True)
    a = ap.parse_args()
    rows = EV.load_rows(a.stage)
    tasks = {t['task']: t for t in json.load(open(os.path.join(WP6, 'tasks', 'tasks.json')))['tasks']}
    audits = {r['run_id']: json.load(open(os.path.join(RUNS, r['run_id'], 'audit.json'))) for r in rows}
    mc_path = os.path.join(RES, f'method_coding_stage{a.stage}.json')
    mc = json.load(open(mc_path)) if os.path.exists(mc_path) else None

    per_task = []
    for t in sorted({r['task'] for r in rows}):
        rr = [r for r in rows if r['task'] == t]
        rel = [abs(r['rel_error']) for r in rr if r['rel_error'] is not None]
        per_task.append(dict(task=t, state_id=tasks[t]['state_id'], task_class=rr[0]['task_class'], G_ref=tasks[t]['reference_G_DD_S_per_cm'],
                             runs=[r['run_id'] for r in rr], outcomes=[r['outcome'] for r in rr],
                             max_abs_rel_error=max(rel) if rel else None, wall_s=[round(r['wall_s']) for r in rr]))
    rel_all = sorted(abs(r['rel_error']) for r in rows if r['rel_error'] is not None)
    rel_summary = dict(n=len(rel_all), median=statistics.median(rel_all) if rel_all else None, max=max(rel_all) if rel_all else None,
                       n_below_1e_3=sum(x <= 1e-3 for x in rel_all), n_below_1e_4=sum(x <= 1e-4 for x in rel_all))
    hints = {k: sorted(r for r, au in audits.items() if au['method_hints'].get(k)) for k in HINT_LABELS}
    methods = None
    if mc:
        methods = {k: dict(label=v, runs=sorted(r for r, x in mc['runs'].items() if x['category'] == k and r in audits))
                   for k, v in mc['categories'].items()}
    ref = {r: refusals(r) for r in sorted(audits)}
    mech = collections.Counter()
    tgt = collections.Counter()
    for r, (m, t) in ref.items():
        mech.update(m)
        tgt.update(t)
    by_class = {}
    for c in ('low', 'moderate'):
        cr = [r for r in rows if r['task_class'] == c]
        by_class[c] = dict(n=len(cr), outcomes=dict(collections.Counter(r['outcome'] for r in cr)),
                           median_wall_s=EV._median([r['wall_s'] for r in cr]), median_turns=EV._median([r['num_turns'] for r in cr if r['num_turns']]),
                           cost_usd=round(sum(r['cost_usd'] or 0 for r in cr), 2))
    out = dict(status='exploratory (post hoc, not pre-registered); descriptive; decides no outcome', stage=a.stage, n_runs=len(rows),
               per_task=per_task, reported_relative_error=rel_summary, by_class=by_class, method_hints=hints, final_method=methods,
               method_coding=os.path.relpath(mc_path, os.path.dirname(os.path.dirname(WP6))) if mc else None,
               refused_calls=dict(by_mechanism=dict(mech), guard_targets=dict(tgt),
                                  runs_with_refusals=sorted(r for r, (m, _) in ref.items() if m), per_run={r: m for r, (m, _) in ref.items() if m}))
    with open(os.path.join(RES, f'exploratory_wp6_stage{a.stage}.json'), 'w') as f:
        json.dump(out, f, indent=1)

    w = print
    n = len(rows)
    w(f'### Exploratory description, stage {a.stage} (post hoc; decides no outcome)\n')
    w('**Per task**\n')
    w('| Task | state | class | G_ref (S/cm) | outcomes | max \\|rel. error\\| | wall (s) |')
    w('|---|---|---|---:|---|---:|---|')
    for x in per_task:
        m = 'n/a' if x['max_abs_rel_error'] is None else f"{x['max_abs_rel_error']:.1e}"
        w(f"| {x['task']} | {x['state_id']} | {x['task_class']} | {x['G_ref']:.4e} | {' '.join(x['outcomes'])} | {m} | {', '.join(map(str, x['wall_s']))} |")
    w(f"\nReported values: {rel_summary['n']}; median |rel. error| {rel_summary['median']:.1e}, largest {rel_summary['max']:.1e}; "
      f"{rel_summary['n_below_1e_4']} within 1e-4 and {rel_summary['n_below_1e_3']} within 1e-3.\n")
    w('**By task class**\n')
    w('| Class | runs | outcomes | median wall (s) | median turns | cost (USD) |')
    w('|---|---:|---|---:|---:|---:|')
    for c, x in by_class.items():
        w(f"| {c} | {x['n']} | {x['outcomes']} | {x['median_wall_s']:.0f} | {x['median_turns']:.0f} | {x['cost_usd']:.2f} |")
    w('\n**Paths the agents\' code touched** (runs whose tool inputs contain the call; mechanical count)\n')
    w('| Path | runs |')
    w('|---|---:|')
    for k, lab in HINT_LABELS.items():
        w(f'| {lab} | {len(hints[k])}/{n} |')
    if methods:
        w('\n**Final method** (analyst coding of each run\'s reported method)\n')
        w('| Category | runs | ids |')
        w('|---|---:|---|')
        for k, x in methods.items():
            w(f"| {k}: {x['label']} | {len(x['runs'])} | {', '.join(x['runs'])} |")
    w('\n**Refused tool calls** (by the mechanism that refused them)\n')
    w('| Mechanism | calls |')
    w('|---|---:|')
    for k, v in sorted(mech.items(), key=lambda kv: -kv[1]):
        w(f'| {k} | {v} |')
    w(f"\nGuard refusals by target: {dict(tgt)}. Runs with any refusal: {', '.join(out['refused_calls']['runs_with_refusals'])}.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
