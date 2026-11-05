"""WP6 harness tests: workspace guard, transcript audit, outcome classification, exact test, prompts, schedule, kit.

They check that the study's instruments behave as pre-registered on synthetic inputs. They are not results.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WP6 = os.path.join(ROOT, 'benchmark', 'wp6')
sys.path.insert(0, os.path.join(WP6, 'harness'))
sys.path.insert(0, os.path.join(WP6, 'agent'))

import audit as AU  # noqa: E402
import evaluate_wp6 as EV  # noqa: E402
import prompts as P  # noqa: E402

GUARD = os.path.join(WP6, 'agent', 'guard.py')
WS = '/home/claude/wp6_runs/W6-999'


def guard(tool, inp):
    d = json.dumps(dict(tool_name=tool, tool_input=inp, cwd=WS))
    return subprocess.run([sys.executable, GUARD], input=d, capture_output=True, text=True).returncode


@pytest.mark.parametrize('tool,inp,code', [
    ('Bash', {'command': './py run.py --vd 1e-4'}, 0),
    ('Bash', {'command': 'ls -la docs/examples'}, 0),
    ('Read', {'file_path': 'docs/devsim.pdf'}, 0),
    ('Read', {'file_path': '/home/claude/venv-wp1/lib/python3.12/site-packages/devsim/__init__.py'}, 0),
    ('Write', {'file_path': f'{WS}/g.py', 'content': 'm = ~np.isnan(x)\n'}, 0),
    ('Bash', {'command': 'ls /home/claude'}, 2),
    ('Bash', {'command': 'cat /home/claude/mosfet-tcad-research/benchmark/README.md'}, 2),
    ('Bash', {'command': 'cd .. && ls'}, 2),
    ('Bash', {'command': 'find / -name "*.json"'}, 2),
    ('Bash', {'command': 'ls ~/.claude'}, 2),
    ('Bash', {'command': 'curl https://example.org'}, 2),
    ('Bash', {'command': 'cat /proc/1234/environ'}, 2),
    ('Read', {'file_path': '/root/.claude/projects/x.jsonl'}, 2),
    ('Grep', {'pattern': 'G_DD', 'path': '/home/claude'}, 2),
    ('Glob', {'pattern': '**/*.json', 'path': '/'}, 2),
    ('Write', {'file_path': f'{WS}/a.py', 'content': 'open("/home/claude/mosfet-tcad-research/x")'}, 2),
])
def test_guard(tool, inp, code):
    assert guard(tool, inp) == code


def _transcript(tmp_path, events, run_id='W6-999', arm='A'):
    rd = tmp_path / 'runs' / run_id
    rd.mkdir(parents=True)
    tr = rd / 'transcript.jsonl.gz'
    with gzip.open(tr, 'wt') as f:
        for e in events:
            f.write(json.dumps(e) + '\n')
    run = dict(run_id=run_id, task='T01', arm=arm, repeat=1, stage=1, probe=False, exit=0, timed_out=False, wall_s=10.0,
               transcript=str(tr))
    (rd / 'run.json').write_text(json.dumps(run))
    return str(rd)


def _ev(final, uses=(), model='claude-opus-5-5', tools=('Bash', 'Read', 'Write', 'Edit', 'Glob', 'Grep')):
    ev = [dict(type='system', subtype='init', tools=list(tools), mcp_servers=[], model=model, permissionMode='auto', cwd=WS)]
    for i, (name, inp, out) in enumerate(uses):
        ev.append(dict(type='assistant', message=dict(model=model, content=[dict(type='tool_use', id=f'u{i}', name=name, input=inp)])))
        ev.append(dict(type='user', message=dict(content=[dict(type='tool_result', tool_use_id=f'u{i}', content=out)])))
    ev.append(dict(type='assistant', message=dict(model=model, content=[dict(type='text', text=final)])))
    ev.append(dict(type='result', subtype='success', result=final, num_turns=3, total_cost_usd=0.5, is_error=False))
    return ev


FINAL = 'Done.\n```json\n{"G_DD_S_per_cm": 3.5e-12, "status": "success", "method": "x", "verification": "none"}\n```'


def test_audit_clean_run_and_final_block(tmp_path, monkeypatch):
    monkeypatch.setattr(AU, 'ROOT', '/')
    rd = _transcript(tmp_path, _ev(FINAL, uses=[('Bash', {'command': './py a.py'}, '3.5e-12'),
                                                ('Bash', {'command': 'ls /home/claude'}, 'Blocked by workspace policy: ...')]))
    a = AU.audit(rd)
    assert a['G_reported'] == 3.5e-12 and a['status'] == 'success'
    assert not a['contaminated'] and not a['harness_failure'] and len(a['blocked_attempts']) == 1


def test_audit_flags_unblocked_access_leak_and_setup(tmp_path, monkeypatch):
    monkeypatch.setattr(AU, 'ROOT', '/')
    rd = _transcript(tmp_path, _ev(FINAL, uses=[('Bash', {'command': 'cat /home/claude/mosfet-tcad-research/x'}, 'secret'),
                                                ('Bash', {'command': './py a.py'}, '"same_mesh_references": {...}')],
                                   tools=('Bash', 'Read', 'Write', 'Edit', 'Glob', 'Grep', 'WebSearch')))
    a = AU.audit(rd)
    assert a['flags']['F1_forbidden_access'] and a['flags']['F2_leaked_material'] and a['flags']['F4_setup']
    assert a['contaminated'] and a['harness_failure']


def test_audit_other_workspace_is_forbidden(tmp_path, monkeypatch):
    monkeypatch.setattr(AU, 'ROOT', '/')
    rd = _transcript(tmp_path, _ev(FINAL, uses=[('Read', {'file_path': '/home/claude/wp6_runs/W6-001/run.py'}, 'x')]))
    assert AU.audit(rd)['contaminated']


@pytest.mark.parametrize('G,status,outcome', [(3.5e-12, 'success', 'CS'), (3.6e-12, 'success', 'SE'), (-3.5e-12, 'success', 'SE'),
                                              (3.6e-12, 'uncertain', 'FE'), (3.5e-12, 'failed', 'FC'), (None, None, 'NR')])
def test_classification(G, status, outcome):
    au = dict(harness_failure=False, contaminated=False, flags={}, G_reported=G, status=status,
              reported=None if G is None else {'G_DD_S_per_cm': G}, checker_calls=[], num_turns=3, n_tool_calls=2,
              total_cost_usd=0.1, result_subtype='success')
    run = dict(run_id='W6-001', task='T01', arm='A', repeat=1, stage=1, wall_s=10.0, timed_out=False)
    assert EV.classify(run, au, 3.5e-12)['outcome'] == outcome


def test_exact_test_matches_brute_force():
    r = EV.selftest()
    assert abs(r['exact_p'] - r['brute_force_p']) < 1e-12


def test_prompts_differ_only_by_checker_paragraph():
    a, b = P.task_prompt(WS, 'A'), P.task_prompt(WS, 'B')
    assert P.CHECKER_PARAGRAPH in b and b.replace(P.CHECKER_PARAGRAPH, '') == a
    assert './check' not in a and './check' in b and 'within 1%' in a


def test_schedule_is_deterministic_and_balanced():
    s = json.load(open(os.path.join(WP6, 'harness', 'schedule.json')))
    runs = [r for r in s['runs'] if r['stage'] in (1, 2)]
    assert len(runs) == 66
    for stage, arm in ((1, 'A'), (2, 'B')):
        rr = [r for r in runs if r['stage'] == stage]
        assert {r['arm'] for r in rr} == {arm}
        assert sorted((r['task'], r['repeat']) for r in rr) == sorted((f'T{k:02d}', j) for k in range(1, 12) for j in (1, 2, 3))


def test_task_files_carry_no_reference_or_identity():
    bad = ('reference_G', 'state_id', 'AR-', 'HO-', 'benchmark', 'network', 'case', 'tolerance', 'G_DD', 'SMOKE')
    tasks = {t['task']: t for t in json.load(open(os.path.join(WP6, 'tasks', 'tasks.json')))['tasks']}
    for t, row in tasks.items():
        ref = f"{row['reference_G_DD_S_per_cm']:.6e}".split('e')[0]
        for f in ('device.json', 'mesh.json'):
            txt = open(os.path.join(WP6, 'tasks', t, f)).read()
            assert not any(b in txt for b in bad), (t, f, [b for b in bad if b in txt])
            assert ref not in txt, (t, f, ref)
