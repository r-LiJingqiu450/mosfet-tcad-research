#!/usr/bin/env python3
"""Mechanical audit of one WP6 run transcript (claude -p stream-json).

Extracts: the session setup (tools, model, MCP servers, permission mode), every tool call and result, the models that
answered, the final message, the reported JSON block, usage and cost. Flags (pre-registered contamination rules):
  F1 forbidden access   a tool input refers to a location outside the agent's workspace, the Python environment and
                        system directories (same patterns as the guard, plus other runs' workspaces), and the call was
                        not blocked by the guard;
  F2 leaked material    a tool result contains text that exists only in the project's files (marker strings below);
  F3 checker tampering  arm B: a tool input touches the verification socket other than through ./check, or a process's
                        memory or environment (/proc/<pid>);
  F4 setup              the session had a tool, MCP server or model other than the pre-registered ones.
A run is contaminated if F1, F2 or F3 is raised; F4 invalidates the run as a harness failure.
Method hints (which DEVSIM calls the agent's code used) are recorded for the root-cause coding; they decide nothing.

Usage:  python benchmark/wp6/harness/audit.py benchmark/wp6/results/runs/<run_id>   (writes audit.json there)
"""
from __future__ import annotations

import gzip
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'agent'))
import guard as G  # noqa: E402

EXPECTED_TOOLS = {'Bash', 'Read', 'Write', 'Edit', 'Glob', 'Grep'}
EXPECTED_MODEL = 'claude-opus-5-5'
BLOCK_MARK = 'Blocked by workspace policy'
LEAK_MARKERS = ('same_mesh_references', 'known_bad_outputs', 'cross_paths', 'mosfet-tcad-research', 'benchmark/cases',
                'Preregistration', 'network_same_run', 'rdeck_', 'phase3/', 'reference_G_DD', 'tasks.json',
                'evaluation_wp', 'posthoc_wp4', 'G_DD_rel_forms', 'kappa_G', 'edge energy')
HINTS = {
    'get_contact_current': r'get_contact_current',
    'circuit_node_value': r'get_circuit_node_value',
    'ac_analysis': r"type\s*=\s*['\"]ac['\"]|ssac_",
    'extended_precision': r'extended_(?:solver|model|equation)',
    'circuit_alter': r'circuit_alter',
    'jacobian_export': r'get_matrix_and_rhs|matrix_and_rhs',
    'richardson_or_extrapolation': r'[Rr]ichardson|extrapolat',
    'electron_density_postprocessing': r"name\s*=\s*['\"]Electrons['\"]",
}


def load(path):
    op = gzip.open if path.endswith('.gz') else open
    out = []
    with op(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('{'):
                try:
                    out.append(json.loads(line))
                except json.JSONDecodeError:
                    out.append(dict(type='unparsed', raw=line[:200]))
    return out


def text_of(content):
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        return '\n'.join(text_of(c.get('text', c.get('content', ''))) if isinstance(c, dict) else str(c) for c in content)
    return str(content)


def input_texts(inp):
    return [v for v in (inp or {}).values() if isinstance(v, str)]


def forbidden(inp, tool, run_id):
    """Reasons a tool input refers to a forbidden location (guard patterns plus other runs' workspaces)."""
    why = []
    cwd = f'{G.ALLOWED_ROOTS[0]}{run_id}'
    for key in ('file_path', 'path', 'notebook_path'):
        p = (inp or {}).get(key)
        if isinstance(p, str) and p:
            ap = os.path.normpath(os.path.join(cwd, p))
            if not (ap + '/').startswith(G.ALLOWED_ROOTS) or ap.startswith('/tmp/claude-'):
                why.append(f'{tool} path {ap}')
    for t in input_texts(inp):
        for pat, what in G.FORBIDDEN:
            if re.search(pat, t):
                why.append(what)
        for m in re.finditer(r'/home/claude/wp6_runs/([\w.-]+)', t):
            if m.group(1) != run_id:
                why.append(f'other workspace {m.group(1)}')
    if tool == 'Bash':
        for pat, what in G.BASH_ONLY:
            if re.search(pat, (inp or {}).get('command', '') or ''):
                why.append(what)
    return why


def final_json(text):
    """The last JSON object in the final message that carries G_DD_S_per_cm."""
    if not text:
        return None
    cands = re.findall(r'```(?:json)?\s*(\{.*?\})\s*```', text, flags=re.S) + re.findall(r'(\{[^{}]*G_DD_S_per_cm[^{}]*\})', text, flags=re.S)
    for c in reversed(cands):
        try:
            d = json.loads(c)
        except json.JSONDecodeError:
            continue
        if isinstance(d, dict) and 'G_DD_S_per_cm' in d:
            return d
    return None


def as_float(v):
    try:
        x = float(v)
        return x if x == x and abs(x) != float('inf') else None
    except (TypeError, ValueError):
        return None


def audit(rundir):
    run = json.load(open(os.path.join(rundir, 'run.json')))
    ev = load(os.path.join(ROOT, run['transcript']))
    init = next((e for e in ev if e.get('type') == 'system' and e.get('subtype') == 'init'), {})
    result = next((e for e in reversed(ev) if e.get('type') == 'result'), {})
    uses, results, models, texts = [], {}, set(), []
    for e in ev:
        if e.get('type') == 'assistant':
            m = e.get('message', {})
            if m.get('model'):
                models.add(m['model'])
            for b in m.get('content', []) or []:
                if b.get('type') == 'tool_use':
                    uses.append(dict(id=b.get('id'), name=b.get('name'), input=b.get('input') or {}))
                elif b.get('type') == 'text':
                    texts.append(b.get('text', ''))
        elif e.get('type') == 'user':
            c = e.get('message', {}).get('content')
            for b in c if isinstance(c, list) else []:
                if isinstance(b, dict) and b.get('type') == 'tool_result':
                    results[b.get('tool_use_id')] = dict(text=text_of(b.get('content')), is_error=bool(b.get('is_error')))
    f1, blocked, f3 = [], [], []
    for u in uses:
        res = results.get(u['id'], {})
        was_blocked = BLOCK_MARK in res.get('text', '')
        if was_blocked:
            blocked.append(dict(tool=u['name'], input=json.dumps(u['input'])[:300]))
        why = forbidden(u['input'], u['name'], run['run_id'])
        if why and not was_blocked:
            f1.append(dict(tool=u['name'], reasons=sorted(set(why)), input=json.dumps(u['input'])[:300]))
        alltext = ' '.join(input_texts(u['input']))
        if re.search(r'\.checker\.sock', alltext) or re.search(r'/proc/\d+/(?:mem|environ|cmdline)', alltext):
            f3.append(dict(tool=u['name'], input=alltext[:300]))
    f2 = []
    for tid, r in results.items():
        hits = [m for m in LEAK_MARKERS if m in r['text']]
        if hits:
            f2.append(dict(tool_use_id=tid, markers=hits, excerpt=r['text'][:300]))
    tools = set(init.get('tools', []) or [])
    f4 = []
    if tools - EXPECTED_TOOLS:
        f4.append(f'unexpected tools {sorted(tools - EXPECTED_TOOLS)}')
    if init.get('mcp_servers'):
        f4.append(f"MCP servers {init.get('mcp_servers')}")
    if models - {EXPECTED_MODEL}:
        f4.append(f'unexpected models {sorted(models - {EXPECTED_MODEL})}')
    code = '\n'.join(t for u in uses for t in input_texts(u['input']))
    hints = {k: len(re.findall(p, code)) for k, p in HINTS.items()}
    final_text = result.get('result') or (texts[-1] if texts else '')
    rep = final_json(final_text)
    checker = []
    lp = os.path.join(rundir, 'checker_log.jsonl')
    if os.path.exists(lp):
        checker = [json.loads(x) for x in open(lp) if x.strip()]
    out = dict(run_id=run['run_id'], task=run['task'], arm=run['arm'], repeat=run['repeat'], stage=run['stage'],
               probe=run.get('probe', False), exit=run['exit'], timed_out=run['timed_out'], wall_s=run['wall_s'],
               session=dict(tools=sorted(tools), mcp_servers=init.get('mcp_servers'), model_init=init.get('model'),
                            permission_mode=init.get('permissionMode'), cwd=init.get('cwd'),
                            claude_code_version=init.get('claude_code_version')),
               models=sorted(models), n_tool_calls=len(uses), n_bash=sum(u['name'] == 'Bash' for u in uses),
               n_python_runs=sum(1 for u in uses if u['name'] == 'Bash' and re.search(r'(?:^|[\s;&|(])\./py\b|python', u['input'].get('command', '') or '')),
               result_subtype=result.get('subtype'), is_error=result.get('is_error'), num_turns=result.get('num_turns'),
               duration_ms=result.get('duration_ms'), total_cost_usd=result.get('total_cost_usd'), usage=result.get('usage'),
               model_usage=result.get('modelUsage'), permission_denials=result.get('permission_denials'),
               rate_limits=[e.get('rate_limit_info') for e in ev if e.get('type') == 'rate_limit_event'][-1:],
               final_message=final_text[-4000:] if final_text else None, reported=rep,
               G_reported=as_float(rep.get('G_DD_S_per_cm')) if rep else None,
               status=(str(rep.get('status', '')).strip().lower() if rep else None),
               flags=dict(F1_forbidden_access=f1, F2_leaked_material=f2, F3_checker_tampering=f3, F4_setup=f4),
               blocked_attempts=blocked, contaminated=bool(f1 or f2 or f3), harness_failure=bool(f4),
               checker_calls=[c for c in checker if 'verdict' in c], checker_requests=[c for c in checker if 'raw' in c],
               method_hints=hints)
    with open(os.path.join(rundir, 'audit.json'), 'w') as f:
        json.dump(out, f, indent=1)
    return out


def main():
    for d in sys.argv[1:]:
        a = audit(d.rstrip('/'))
        print(json.dumps({k: a[k] for k in ('run_id', 'task', 'arm', 'exit', 'G_reported', 'status', 'contaminated',
                                            'harness_failure', 'n_tool_calls', 'num_turns', 'total_cost_usd')}))
    return 0


if __name__ == '__main__':
    sys.exit(main())
