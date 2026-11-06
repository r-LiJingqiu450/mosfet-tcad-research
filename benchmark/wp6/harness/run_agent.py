#!/usr/bin/env python3
"""Run one WP6 agent session (one run of the schedule) and archive everything it produced.

Steps:
  1. build a fresh workspace /home/claude/wp6_runs/<run_id> from the kit and the task's data (no reference, no case id);
  2. arm B only: start the verification service (checker_daemon.py), reference passed on its standard input;
  3. run the agent under test as a fresh headless Claude Code session (`claude -p`), with
     model claude-opus-5-5, effort high, tools Bash/Read/Write/Edit/Glob/Grep only, no MCP servers (also no Claude in
     Chrome), no skills or slash commands, no settings files, permission mode auto (Claude Code's action classifier
     reviews every tool call), the workspace guard as PreToolUse hook, the fixed system prompt (prompts.SYSTEM_PROMPT),
     at most MAX_TURNS turns and WALL_S seconds of wall time, full stream-json transcript;
  4. archive transcript, stderr, checker log and the workspace files the agent created into
     benchmark/wp6/results/runs/<run_id>/, append one line to benchmark/wp6/results/ledger.jsonl, delete the workspace.

Usage (repository root):  python benchmark/wp6/harness/run_agent.py --run-id W6-001 [--probe "text"]
The run's task, arm and repeat come from benchmark/wp6/harness/schedule.json (smoke runs: S1, S2). --probe replaces the
task prompt (harness checks only; such runs are labelled probe and are never part of the study).
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prompts as P  # noqa: E402

WP6 = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(WP6))
KIT = os.path.join(WP6, 'kit')
TASKS = os.path.join(WP6, 'tasks')
RESULTS = os.path.join(WP6, 'results')
WS_ROOT = '/home/claude/wp6_runs'
GUARD = '/opt/wp6/guard.py'
MODEL, EFFORT, MAX_TURNS, WALL_S, PERMISSION_MODE = 'claude-opus-5-5', 'high', 150, 2700, 'auto'
TOOLS = 'Bash,Read,Write,Edit,Glob,Grep'
SETTINGS = {'hooks': {'PreToolUse': [{'matcher': 'Bash|Read|Write|Edit|Glob|Grep',
                                      'hooks': [{'type': 'command', 'command': f'python3 {GUARD}'}]}]}}
KIT_FILES = ('mosfet_deck.py', 'py')
TASK_FILES = ('device.json', 'mesh.json', 'zero_bias_state.npz')
# environment variables of this (parent) session that the agent under test does not get: messaging and session
# identifiers of the parent, account identifiers, and files of the parent session
DROP_ENV_PREFIXES = ('CLAUDE_CODE_MESSAGING', 'CLAUDE_CODE_SESSION_ID', 'CLAUDE_CODE_REMOTE_SESSION_ID', 'CLAUDE_PROJECT_',
                     'CLAUDE_CODE_USER_EMAIL', 'CLAUDE_CODE_ACCOUNT_UUID', 'CLAUDE_CODE_ORGANIZATION_UUID',
                     'CLAUDE_CODE_DIAGNOSTICS_FILE', 'CLAUDE_CODE_DEBUG', 'CLAUDE_AFTER_LAST_COMPACT', 'CLAUDE_PID',
                     'CLAUDE_CODE_CHILD_SESSION', 'CLAUDE_CODE_ADDITIONAL_DIRECTORIES_CLAUDE_MD', 'CLAUDE_ADDITIONAL_DIRECTORIES',
                     'CLAUDE_CODE_ARTIFACT', 'CLAUDECODE', 'CLAUDE_CHROME', 'CLAUDE_CODE_CHROME',
                     'CLAUDE_CODE_INCLUDE_PARTIAL_MESSAGES', 'CLAUDE_CODE_ENABLE_PROMPT_SUGGESTION')


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def sha_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def load_schedule():
    return json.load(open(os.path.join(HERE, 'schedule.json')))


def task_row(task):
    for t in json.load(open(os.path.join(TASKS, 'tasks.json')))['tasks']:
        if t['task'] == task:
            return t
    raise SystemExit(f'unknown task {task}')


def make_workspace(run_id, task, arm):
    ws = os.path.join(WS_ROOT, run_id)
    if os.path.exists(ws):
        raise SystemExit(f'workspace exists: {ws}')
    os.makedirs(ws)
    for f in KIT_FILES:
        shutil.copy2(os.path.join(KIT, f), ws)
    for f in TASK_FILES:
        shutil.copy2(os.path.join(TASKS, task, f), ws)
    shutil.copytree(os.path.join(KIT, 'docs'), os.path.join(ws, 'docs'))
    if arm == 'B':
        shutil.copy2(os.path.join(KIT, 'check'), ws)
    with open(os.path.join(ws, 'TASK.md'), 'w') as f:
        f.write(P.task_prompt(ws, arm) + '\n')
    manifest = {}
    for d, _, fs in os.walk(ws):
        for f in fs:
            p = os.path.join(d, f)
            manifest[os.path.relpath(p, ws)] = sha_file(p)
    return ws, manifest


def child_env():
    return {k: v for k, v in os.environ.items() if not k.startswith(DROP_ENV_PREFIXES)}


def kill_strays(ws):
    """Kill processes still running in the workspace after the session ended (e.g. a simulation left behind)."""
    killed = []
    for pid in os.listdir('/proc'):
        if not pid.isdigit() or int(pid) == os.getpid():
            continue
        try:
            cwd = os.readlink(f'/proc/{pid}/cwd')
        except OSError:
            continue
        if cwd == ws or cwd.startswith(ws + '/'):
            try:
                cmd = open(f'/proc/{pid}/cmdline', 'rb').read().replace(b'\0', b' ').decode(errors='replace')[:200]
                os.kill(int(pid), signal.SIGKILL)
                killed.append(dict(pid=int(pid), cmd=cmd))
            except OSError:
                pass
    return killed


def snapshot(ws, dst, manifest):
    """Copy every file the agent created or changed (files under 20 MB), skip the unchanged kit files."""
    out = []
    for d, _, fs in os.walk(ws):
        for f in fs:
            p = os.path.join(d, f)
            rel = os.path.relpath(p, ws)
            if rel == '.checker.sock' or (rel in manifest and os.path.exists(p) and sha_file(p) == manifest[rel]):
                continue
            size = os.path.getsize(p) if os.path.isfile(p) else 0
            rec = dict(path=rel, size=size)
            if os.path.isfile(p) and size <= 20e6:
                os.makedirs(os.path.dirname(os.path.join(dst, rel)), exist_ok=True)
                shutil.copy2(p, os.path.join(dst, rel))
                rec['sha256'] = sha_file(p)
            else:
                rec['skipped'] = 'larger than 20 MB or not a regular file'
            out.append(rec)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--probe', help='replace the task prompt (harness check; not part of the study)')
    a = ap.parse_args()
    sched = {r['run_id']: r for r in load_schedule()['runs']}
    r = sched[a.run_id] if not a.probe else dict(run_id=a.run_id, task='SMOKE', arm='A', repeat=0, stage='probe')
    t = task_row(r['task'])
    rundir = os.path.join(RESULTS, 'runs', a.run_id)
    if os.path.exists(os.path.join(rundir, 'run.json')):
        raise SystemExit(f'{a.run_id} already has a record')
    os.makedirs(rundir, exist_ok=True)
    ws, manifest = make_workspace(a.run_id, r['task'], r['arm'])
    prompt = a.probe if a.probe else P.task_prompt(ws, r['arm'])
    rec = dict(run_id=a.run_id, task=r['task'], arm=r['arm'], repeat=r['repeat'], stage=r['stage'], probe=bool(a.probe),
               start_utc=now(), workspace=ws, model=MODEL, effort=EFFORT, max_turns=MAX_TURNS, wall_limit_s=WALL_S,
               permission_mode=PERMISSION_MODE,
               tools=TOOLS, settings=SETTINGS, guard_sha256=sha_file(GUARD),
               system_prompt_sha256=P.sha(P.SYSTEM_PROMPT), task_prompt_sha256=P.sha(prompt),
               workspace_manifest=manifest, claude_version=subprocess.run(['claude', '--version'], capture_output=True,
                                                                          text=True).stdout.strip())
    daemon = None
    if r['arm'] == 'B':
        daemon = subprocess.Popen([sys.executable, os.path.join(HERE, 'checker_daemon.py'), '--socket',
                                   os.path.join(ws, '.checker.sock'), '--log', os.path.join(rundir, 'checker_log.jsonl')],
                                  stdin=subprocess.PIPE, text=True)
        daemon.stdin.write(repr(t['reference_G_DD_S_per_cm']) + '\n'); daemon.stdin.close()
        for _ in range(50):
            if os.path.exists(os.path.join(ws, '.checker.sock')):
                break
            time.sleep(0.1)
    cmd = ['timeout', '-k', '60', str(WALL_S), 'claude', '-p', '--model', MODEL, '--effort', EFFORT, '--tools', TOOLS,
           '--strict-mcp-config', '--setting-sources', '', '--settings', json.dumps(SETTINGS),
           '--no-chrome', '--disable-slash-commands', '--system-prompt', P.SYSTEM_PROMPT, '--permission-mode', PERMISSION_MODE,
           '--output-format', 'stream-json', '--verbose', '--no-session-persistence', '--max-turns', str(MAX_TURNS)]
    t0 = time.time()
    tr = os.path.join(rundir, 'transcript.jsonl')
    with open(tr, 'w') as out, open(os.path.join(rundir, 'stderr.log'), 'w') as err:
        p = subprocess.run(cmd, input=prompt, stdout=out, stderr=err, cwd=ws, env=child_env(), text=True)
    rec.update(end_utc=now(), wall_s=round(time.time() - t0, 1), exit=p.returncode, timed_out=p.returncode in (124, 137))
    rec['stray_processes_killed'] = kill_strays(ws)
    if daemon is not None:
        daemon.terminate()
        daemon.wait(timeout=10)
    rec['created_files'] = snapshot(ws, os.path.join(rundir, 'workspace'), manifest)
    shutil.rmtree(ws)
    with open(tr, 'rb') as f, gzip.open(tr + '.gz', 'wb', compresslevel=9) as g:
        shutil.copyfileobj(f, g)
    os.remove(tr)
    rec['transcript'] = os.path.relpath(tr + '.gz', ROOT)
    rec['transcript_sha256'] = sha_file(tr + '.gz')
    with open(os.path.join(rundir, 'run.json'), 'w') as f:
        json.dump(rec, f, indent=1)
    with open(os.path.join(RESULTS, 'ledger.jsonl'), 'a') as f:
        f.write(json.dumps({k: rec[k] for k in ('run_id', 'task', 'arm', 'repeat', 'stage', 'probe', 'start_utc', 'end_utc',
                                                 'wall_s', 'exit', 'timed_out')}) + '\n')
    print(f"{a.run_id} task={r['task']} arm={r['arm']} exit={p.returncode} wall={rec['wall_s']}s", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
