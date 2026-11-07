"""Release integrity: every file a case cites lives at its recorded path with its recorded SHA-256, the benchmark
modules import from the files present, and (inside a release snapshot) SHA256SUMS matches every file.

These tests run both in the research repository and in the release snapshot. They check the package, not physics.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def sha(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def test_every_case_source_is_present_with_its_sha256():
    n = 0
    for p in sorted(glob.glob(os.path.join(ROOT, 'benchmark', 'cases', '*.json'))):
        for s in json.load(open(p))['provenance'].get('sources', []):
            path = os.path.join(ROOT, s['path'])
            assert os.path.exists(path), (os.path.basename(p), s['path'])
            assert sha(path) == s['sha256'], (os.path.basename(p), s['path'])
            n += 1
    assert n > 0


@pytest.mark.parametrize('module', ['benchmark/tools/score.py', 'benchmark/wp2/rdeck_common.py', 'benchmark/wp2/fv_rdeck.py',
                                    'benchmark/wp8/ssmech.py', 'benchmark/wp8/xsens.py'])
def test_benchmark_modules_resolve_their_imports(module):
    """Import (not run) a module in a fresh interpreter, with only the release root on the path."""
    if 'rdeck_common' in module or '/wp8/' in module:
        os.environ.setdefault('DEVSIM_MATH_LIBS', '/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3')
        pytest.importorskip('devsim', reason='DEVSIM not installed')
    code = ('import importlib.util, sys; sys.path[:0] = [sys.argv[1], sys.argv[2]]; '
            'spec = importlib.util.spec_from_file_location("m", sys.argv[3]); m = importlib.util.module_from_spec(spec); '
            'spec.loader.exec_module(m)')
    env = dict(os.environ)
    env.setdefault('DEVSIM_MATH_LIBS', '/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3')
    p = subprocess.run([sys.executable, '-c', code, ROOT, os.path.join(ROOT, os.path.dirname(module)),
                        os.path.join(ROOT, module)], cwd=ROOT, env=env, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]


def test_release_checksums_when_in_a_snapshot():
    sums = os.path.join(ROOT, 'SHA256SUMS')
    if not os.path.exists(sums):
        pytest.skip('not a release snapshot')
    lines = [x.split('  ', 1) for x in open(sums).read().splitlines() if x]
    listed = {path for _, path in lines}
    for h, path in lines:
        assert sha(os.path.join(ROOT, path)) == h, path
    present = {os.path.relpath(os.path.join(dp, f), ROOT) for dp, _, fs in os.walk(ROOT) for f in fs}
    # files made by running the code here (bytecode caches, git metadata) are not release content unless listed
    present = {p for p in present if p != 'SHA256SUMS' and (p in listed or not (
        p.startswith(('.git/', '.pytest_cache/')) or '__pycache__' in p))}
    assert present == listed, (sorted(present - listed)[:10], sorted(listed - present)[:10])
