#!/usr/bin/env python3
"""Install the WP6 workspace guard (PreToolUse hook of the agent under test) to /opt/wp6/guard.py.

Prints its SHA-256; every run record carries it. Usage (repository root):  python benchmark/wp6/harness/install_guard.py
"""
from __future__ import annotations

import hashlib
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), 'agent', 'guard.py')
DST = '/opt/wp6/guard.py'


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def main():
    os.makedirs(os.path.dirname(DST), exist_ok=True)
    shutil.copyfile(SRC, DST)
    os.chmod(DST, 0o755)
    assert sha(DST) == sha(SRC)
    print(f'{DST}  sha256={sha(DST)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
