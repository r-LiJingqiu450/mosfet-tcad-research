"""C-10 probe: phase3/c7/endprobe.py with a cancellation-free evaluation of the zero-bias network conductance.

Numerical issue E10-1 (found in C-10):
- wprobe.solve_net imposes phi = 0 on the source contact and phi = 1 on the drain contact.
- It returns G as the DRAIN current, I_d = sum over drain nodes of c_ij (1 - phi_j).
- In deep subthreshold the edge conductances span more than 30 orders of magnitude, and near the drain
  1 - phi_j is far below the double-precision spacing of numbers near 1. The drain-side sum then loses all accuracy:
  it violated the Rayleigh bounds by up to 120x.
- The SOURCE current is a sum of positive terms c_ij phi_j with phi_j close to 0, which double precision represents
  with full relative accuracy.

What this probe changes:
- G is the source-side current of the forward solve (phi_s = 0, phi_d = 1).
- As an independent check, the swapped problem (phi_s = 1, phi_d = 0) is solved and its drain-side current, now the
  near-zero side, is recorded. By reciprocity the two must agree.

Physics, mesh and every other output are unchanged from endprobe.
Diagnostics go to <out>/netcheck_<state>.json; the state JSON is written by endprobe as usual.
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np  # noqa: E402

from phase3.c6 import wprobe as WP  # noqa: E402
from phase3.c7 import endprobe as E7  # noqa: E402

_SOLVE0 = WP.solve_net
NETDIAG = []


def solve_net_robust(N, n0i, n1i, c, src, drn):
    G_d, I_s, phi, res = _SOLVE0(N, n0i, n1i, c, src, drn)          # forward: phi_s = 0, phi_d = 1
    G_d2, I_s2, phi2, res2 = _SOLVE0(N, n0i, n1i, c, drn, src)      # swapped: phi_d = 0, phi_s = 1
    # in the swapped solve the near-zero side is the drain (passed as 'src'); its current is I_s2
    G_fwd_src, G_swp_drn = -I_s, -I_s2
    NETDIAG.append(dict(G_drain_side_forward=G_d, G_source_side_forward=G_fwd_src, G_near_zero_side_swapped=G_swp_drn,
                        rel_src_vs_swapped=G_fwd_src / G_swp_drn - 1, rel_drainside_error=G_d / G_fwd_src - 1, res=res, res_swapped=res2))
    return G_fwd_src, I_s, phi, res


def main():
    WP.solve_net = solve_net_robust
    out = sys.argv[sys.argv.index('--out') + 1]
    E7.main()
    a = E7.argparse  # noqa: F841  (argparse already consumed by E7.main)
    tag = '_'.join(x.lstrip('-') for x in sys.argv[1:] if x != out and x != '--out')
    json.dump(dict(argv=sys.argv[1:], networks=['const_mu', 'lombardi_mu'][:len(NETDIAG)], diag=NETDIAG),
              open(os.path.join(out, f'netcheck_{tag}.json'), 'w'), indent=1)
    for d in NETDIAG:
        print(f"G(source side) {d['G_source_side_forward']:.6e}  swapped {d['G_near_zero_side_swapped']:.6e}  "
              f"rel {d['rel_src_vs_swapped']:.2e}  drain-side error {d['rel_drainside_error']:.2e}")


if __name__ == '__main__':
    main()
