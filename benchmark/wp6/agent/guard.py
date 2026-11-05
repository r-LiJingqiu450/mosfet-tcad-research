#!/usr/bin/env python3
"""WP6 workspace guard: PreToolUse hook of the agent under test (tcad-task-agent).

Reads the hook input (JSON on stdin) and blocks a tool call (exit code 2, reason on stderr) when it refers to a location
outside the agent's workspace, the Python environment and system directories. It is a first barrier only; the binding
control is the transcript audit (benchmark/wp6/harness/audit.py), which also sees what the guard cannot (for example
paths built at run time). Installed to /opt/wp6/guard.py by benchmark/wp6/harness/install_agent.py.
"""
from __future__ import annotations

import json
import os
import re
import sys

ALLOWED_ROOTS = ('/home/claude/wp6_runs/', '/home/claude/venv-wp1/', '/usr/', '/lib/', '/lib64/', '/etc/', '/tmp/', '/dev/null')
# text patterns refused anywhere in a tool input (commands, file contents, search patterns)
FORBIDDEN = [
    (r'/home(?!/claude/(?:wp6_runs|venv-wp1)(?:/|\b))', 'locations under /home other than the workspace and /home/claude/venv-wp1'),
    (r'/tmp/claude-', 'that temporary directory'),
    (r'(?<![\w.~-])/root(?:/|\b)', '/root'),
    (r'/var/lib/wp6|/opt/wp6', 'harness files'),
    (r'/proc/(?!self\b|cpuinfo|meminfo)', 'other processes'),
    (r'expanduser|Path\.home\(\)', 'the home directory'),
]
BASH_ONLY = [
    (r'(?:^|[\s=:;|&(])~(?:/|\s|$|;)', 'the home directory'),
    (r'\$\{?HOME\b', 'the home directory'),
    (r'(?:^|[\s;|&(])(?:find|du|tree|locate|updatedb)\s+/(?:\s|$)', 'searching the whole file system'),
    (r'(?:^|[\s;|&(])(?:ls|find|du|tree|grep\s+-\w*r\w*)\s+(?:-\S+\s+)*(?:/|/home|/home/claude)/?(?:\s|$|;)', 'listing other locations'),
    (r'(?:\.\./){2,}|(?:^|[\s;|&(])cd\s+\.\.(?:\s|$|;|/)', 'leaving the workspace'),
    (r'\b(?:curl|wget|scp|ssh|rsync)\b|\bgit\s+(?:clone|fetch|pull)\b|https?://|\bpip3?\s+install\b|\buv\s+pip\b', 'network access or installing packages'),
]


def refuse(why):
    sys.stderr.write(f'Blocked by workspace policy: {why} is not allowed. Work inside your workspace; you may read the Python '
                     'environment at /home/claude/venv-wp1 and standard system files.\n')
    sys.exit(2)


def main():
    try:
        d = json.load(sys.stdin)
    except Exception:  # noqa: BLE001 - malformed input: do not block
        return 0
    tool = d.get('tool_name', '')
    inp = d.get('tool_input') or {}
    cwd = d.get('cwd') or '/'
    for key in ('file_path', 'path', 'notebook_path'):
        p = inp.get(key)
        if isinstance(p, str) and p:
            ap = os.path.normpath(os.path.join(cwd, os.path.expanduser(p)))
            if not (ap + '/').startswith(ALLOWED_ROOTS) or ap.startswith('/tmp/claude-'):
                refuse(f'{tool} on {ap}')
    texts = [v for v in inp.values() if isinstance(v, str)]
    for t in texts:
        for pat, what in FORBIDDEN:
            if re.search(pat, t):
                refuse(f'access to {what}')
    if tool == 'Bash':
        cmd = inp.get('command', '') or ''
        for pat, what in BASH_ONLY:
            if re.search(pat, cmd):
                refuse(what)
    return 0


if __name__ == '__main__':
    sys.exit(main())
