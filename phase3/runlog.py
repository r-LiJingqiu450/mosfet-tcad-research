"""Run-record helper (Research_05 v2.1 rule 3): config, git SHA, versions, timestamps."""
from __future__ import annotations

import datetime as _dt
import json
import os
import platform
import subprocess
import sys


def _git(*args):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        return subprocess.run(['git', '-C', root, *args], capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception as exc:  # pragma: no cover
        return f'unavailable: {exc}'


def versions():
    out = dict(python=sys.version.split()[0], executable=sys.executable,
               platform=platform.platform())
    for mod in ('numpy', 'scipy', 'matplotlib', 'devsim'):
        try:
            m = __import__(mod)
            out[mod] = getattr(m, '__version__', 'unknown')
        except Exception as exc:
            out[mod] = f'not importable: {type(exc).__name__}'
    out['DEVSIM_MATH_LIBS'] = os.environ.get('DEVSIM_MATH_LIBS', '')
    return out


def now():
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec='milliseconds')


def header(config: dict):
    return dict(start_utc=now(), git_sha=_git('rev-parse', 'HEAD'),
                git_branch=_git('rev-parse', '--abbrev-ref', 'HEAD'),
                git_dirty=bool(_git('status', '--porcelain', '--untracked-files=no')),
                argv=sys.argv, versions=versions(), config=config)


def _clean(o):
    if isinstance(o, dict):
        return {(':'.join(map(str, k)) if isinstance(k, tuple) else k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def _default(o):
    try:
        return float(o)
    except Exception:
        return str(o)


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(_clean(obj), fh, indent=2, default=_default)
        fh.write('\n')
