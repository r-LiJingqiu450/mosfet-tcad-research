#!/usr/bin/env python3
"""WP4-B: the WP2 drift-diffusion path P2 (benchmark/wp2/rdeck_dd.py, unchanged) with DEVSIM extended precision.

DEVSIM's documented switches `extended_solver`, `extended_model` and `extended_equation` are set before anything is
built; then rdeck_dd.main() runs with the given arguments (same mesh, same equilibrium procedure, same Newton rule
tau = 1e-12 V with 2 extra iterations, same four biases and the same four current evaluations). In DEVSIM 2.11.2 the
Newton residual, models and equations are then evaluated in float128 while the linear solves still use the double
precision UMFPACK factorization (src/math/ExternalPreconditioner.cc converts to double): mixed-precision Newton.

Output: <out>/P2X_<id>_<body>_lev<level>.json = the P2 record plus a 'precision' block (DEVSIM build info and the three
parameters read back after the run).

Usage (repository root):  python benchmark/wp4/rdeck_dd_ext.py --id ... --na ... --tox ... --L ... --vg ... --level N
                          --body B1 --out DIR
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WP2 = os.path.join(os.path.dirname(HERE), 'wp2')
sys.path.insert(0, WP2)
from rdeck_common import ds  # noqa: E402  (imports devsim; nothing is built yet)

FLAGS = ('extended_solver', 'extended_model', 'extended_equation')
WRAPPER_VERSION = 'wp4-rdeck-dd-ext-2026-11-04'


def main():
    info = ds.get_parameter(name='info')
    if not info.get('extended_precision'):
        raise SystemExit('this DEVSIM build has no extended precision')
    for p in FLAGS:
        ds.set_parameter(name=p, value=True)
    import rdeck_dd  # noqa: E402
    rc = rdeck_dd.main()
    a = rdeck_dd.argparse.Namespace(**{k: v for k, v in zip(sys.argv[1::2], sys.argv[2::2])})
    out_dir = getattr(a, '--out')
    tag = f"{getattr(a, '--id')}_{getattr(a, '--body')}_lev{getattr(a, '--level')}"
    src, dst = os.path.join(out_dir, f'P2_{tag}.json'), os.path.join(out_dir, f'P2X_{tag}.json')
    rec = json.load(open(src))
    rec['precision'] = dict(wrapper=WRAPPER_VERSION, devsim_info={k: (list(v) if isinstance(v, tuple) else v) for k, v in info.items()},
                            parameters_after_run={p: bool(ds.get_parameter(name=p)) for p in FLAGS},
                            note='Newton residual, models and equations in float128; linear solves in double (UMFPACK)')
    with open(dst, 'w', encoding='utf-8') as f:
        json.dump(rec, f, indent=1)
    os.remove(src)
    print(f'wrote {dst}', flush=True)
    return rc


if __name__ == '__main__':
    sys.exit(main())
