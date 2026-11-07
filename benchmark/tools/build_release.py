#!/usr/bin/env python3
"""Build the public release snapshot of the low-current MOSFET drift-diffusion benchmark (WP7).

The snapshot keeps every file at its repository path, so that every SHA-256 a case file records and every import of
the benchmark code resolve inside it. It holds:
  - benchmark/ : every tracked file;
  - the result files outside benchmark/ that case files cite as sources (WP1 exports: results/phase3/...);
  - the project modules that the benchmark code imports or that case files name as producing scripts (CODE below). The
    list was found in WP7 by importing every module of benchmark/wp2, benchmark/wp4 and benchmark/tools and recording the
    repository files loaded, plus a static import pass over the producing scripts. The snapshot tests import them again;
  - paper/ : manuscript source, the scripts that make its figures, tables and numbers, and the built PDF;
  - release/* copied to the snapshot root: README.md, LICENSE, LICENSE-DATA.md, CITATION.cff, .zenodo.json,
    RELEASE_NOTES.md.
Only committed content is used (git ls-files and git show HEAD:<path>), so the snapshot is the state of one commit.

Writes dist/<name>-<version>/ with SHA256SUMS and a reproducible <name>-<version>.tar.gz with its .sha256. Optionally
it commits the snapshot (a tree holding exactly the snapshot files):
  --branch                 as the single commit of a new orphan branch release/v<version>;
  --update-branch NAME     as a new commit on top of the existing branch NAME (a fast-forward: the earlier snapshot
                           commits stay in its history; the ref is moved only if it still points at the commit read).
The commit message records the source branch and commit and the tarball's SHA-256; --trailer adds trailer lines.
The current checkout and index are not touched.

Usage (repository root):
  python benchmark/tools/build_release.py --version 1.0.0 [--branch | --update-branch release/v1.0.0-candidate]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NAME = 'mosfet-lowcurrent-dd-benchmark'
CODE = [
    'phase3/__init__.py', 'phase3/m0.py', 'phase3/meshmap.py', 'phase3/runlog.py',
    'phase3/c1/__init__.py', 'phase3/c1/c1b_jacobian.py', 'phase3/c1/corner_study.py', 'phase3/c1/deck.py',
    'phase3/c1/measure.py', 'phase3/c1/onecol.py',
    'phase3/c3/__init__.py', 'phase3/c3/t1.py',
    'phase3/c4/__init__.py', 'phase3/c4/deck4.py', 'phase3/c4/t1v2.py',
    'phase3/c6/__init__.py', 'phase3/c6/rprobe.py', 'phase3/c6/wprobe.py',
    'phase3/c7/__init__.py', 'phase3/c7/endprobe.py',
    'phase3/c11/__init__.py', 'phase3/c11/fd2d_indep.py',
    'phase3/direction_reset/well_depth_probe.py', 'phase3/next_question/physical_r2_finite_bias.py',
    'semiconductor_sim/common/__init__.py', 'semiconductor_sim/common/constants.py',   # imported as `common` (phase3/__init__.py)
]
TOP = ['README.md', 'LICENSE', 'LICENSE-DATA.md', 'CITATION.cff', '.zenodo.json', 'RELEASE_NOTES.md']


def git(*args, binary=False):
    out = subprocess.run(['git', *args], cwd=ROOT, capture_output=True, check=True)
    return out.stdout if binary else out.stdout.decode()


def tracked(prefix):
    return [p for p in git('ls-files', '-z', prefix).split('\0') if p]


def sha256(b):
    return hashlib.sha256(b).hexdigest()


def case_sources():
    """(path, sha256) of every case-file source outside benchmark/."""
    out = {}
    for p in tracked('benchmark/cases'):
        if not p.endswith('.json'):
            continue
        case = json.loads(git('show', f'HEAD:{p}'))
        for s in case['provenance'].get('sources', []):
            if not s['path'].startswith('benchmark/'):
                if s['path'] in out and out[s['path']] != s['sha256']:
                    raise SystemExit(f"{s['path']}: two case files record different SHA-256")
                out[s['path']] = s['sha256']
    return out


def plan():
    """Map snapshot path -> repository path (all committed)."""
    files = {p: p for p in tracked('benchmark') + tracked('paper')}
    srcs = case_sources()
    all_tracked = set(git('ls-files', '-z').split('\0'))
    for p in list(srcs) + CODE:
        if p not in all_tracked:
            raise SystemExit(f'{p}: needed by the snapshot but not tracked in git')
        files[p] = p
    for t in TOP:
        rp = f'release/{t}'
        if rp not in all_tracked:
            raise SystemExit(f'{rp}: missing (release files must be committed)')
        files[t] = rp
    return files, srcs


def build(version, branch, update_branch=None, trailers=()):
    files, srcs = plan()
    out = os.path.join(ROOT, 'dist', f'{NAME}-{version}')
    if os.path.exists(out):
        shutil.rmtree(out)
    sums = []
    for dst, src in sorted(files.items()):
        b = git('show', f'HEAD:{src}', binary=True)
        p = os.path.join(out, dst)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, 'wb') as f:
            f.write(b)
        if src in srcs and sha256(b) != srcs[src]:
            raise SystemExit(f'{src}: SHA-256 differs from the one its case files record')
        sums.append(f'{sha256(b)}  {dst}')
    with open(os.path.join(out, 'SHA256SUMS'), 'w') as f:
        f.write('\n'.join(sums) + '\n')
    # reproducible tarball: sorted members, fixed mtime/owner, gzip without timestamp
    tar_path = os.path.join(ROOT, 'dist', f'{NAME}-{version}.tar.gz')
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode='w', format=tarfile.PAX_FORMAT) as tar:
        for dst in sorted(list(files) + ['SHA256SUMS']):
            p = os.path.join(out, dst)
            ti = tar.gettarinfo(p, arcname=f'{NAME}-{version}/{dst}')
            ti.mtime, ti.uid, ti.gid, ti.uname, ti.gname, ti.mode = 0, 0, 0, '', '', 0o644
            with open(p, 'rb') as f:
                tar.addfile(ti, f)
    with open(tar_path, 'wb') as f:
        with gzip.GzipFile(filename='', mode='wb', fileobj=f, mtime=0, compresslevel=9) as gz:
            gz.write(raw.getvalue())
    with open(tar_path, 'rb') as f:
        tsha = sha256(f.read())
    with open(tar_path + '.sha256', 'w') as f:
        f.write(f'{tsha}  {os.path.basename(tar_path)}\n')
    commit, ref, parent, unchanged = None, None, None, False
    if branch or update_branch:
        env = dict(os.environ, GIT_INDEX_FILE=os.path.join(ROOT, 'dist', '.release-index'), GIT_WORK_TREE=out)
        if os.path.exists(env['GIT_INDEX_FILE']):
            os.remove(env['GIT_INDEX_FILE'])
        subprocess.run(['git', 'add', '-A', '-f', '.'], cwd=out, env=dict(env, GIT_DIR=os.path.join(ROOT, '.git')), check=True)
        tree = subprocess.run(['git', 'write-tree'], cwd=out, env=dict(env, GIT_DIR=os.path.join(ROOT, '.git')),
                              capture_output=True, text=True, check=True).stdout.strip()
        os.remove(env['GIT_INDEX_FILE'])
        src_commit = git('rev-parse', 'HEAD').strip()
        src_ref = git('rev-parse', '--abbrev-ref', 'HEAD').strip()
        msg = (f'{NAME} v{version}\n\nRelease snapshot built by benchmark/tools/build_release.py from the research repository '
               f'commit {src_commit} (branch {src_ref}).\nTarball {os.path.basename(tar_path)}: SHA-256 {tsha}.\n'
               f'Code: MIT (LICENSE). Data, reports and paper text: CC BY 4.0 (LICENSE-DATA.md).\n')
        if update_branch:
            ref = f'refs/heads/{update_branch}'
            parent = git('rev-parse', '--verify', f'{ref}^{{commit}}').strip()
            if git('rev-parse', f'{parent}^{{tree}}').strip() == tree:
                unchanged = True
            else:
                msg += f'\nUpdates the earlier snapshot {parent}.\n'
        else:
            ref = f'refs/heads/release/v{version}'
        if trailers:
            msg += '\n' + '\n'.join(trailers) + '\n'
        if not unchanged:
            cmd = ['git', 'commit-tree', tree] + (['-p', parent] if parent else []) + ['-m', msg]
            commit = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
            # compare-and-swap: refuse if the branch moved since it was read (new branch: must not exist yet)
            git('update-ref', ref, commit, parent if parent else '0' * 40)
    print(json.dumps(dict(snapshot=os.path.relpath(out, ROOT), files=len(files) + 1, case_sources_outside_benchmark=len(srcs),
                          code_files=len(CODE), tarball=os.path.relpath(tar_path, ROOT), tarball_sha256=tsha,
                          branch=ref[len('refs/heads/'):] if ref else None, branch_commit=commit, parent=parent,
                          tree_unchanged=unchanged,
                          size_mb=round(sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(out) for f in fs) / 1e6, 1)),
                     indent=1))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True)
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--branch', action='store_true', help='also create the orphan branch release/v<version>')
    g.add_argument('--update-branch', metavar='NAME', help='also commit the snapshot on top of the existing branch NAME')
    ap.add_argument('--trailer', action='append', default=[], help='trailer line for the commit message (repeatable)')
    a = ap.parse_args()
    return build(a.version, a.branch, a.update_branch, a.trailer)


if __name__ == '__main__':
    sys.exit(main())
