"""C-9 post-hoc mesh diagnostic (NOT pre-registered; disclosed as such in the C-9 report).

Question: is the L = 0.5 um deviation of the Lombardi-mobility end correction in the 2e18-halo device
(dL_mu(0.5 um) differs from dL_mu(1 um) = dL_mu(2 um) by up to ~170 nm, while the constant-mobility dL agrees
to 1.2 nm) a length effect or a mesh effect?  For L < 0.7 um, probe9 shrinks the graded lateral end segment to
a = 0.35 L, which coarsens the channel mesh in the halo region (h ~ 2.8-6.9 nm at 20-50 nm from the junction,
level 0, versus 1.6-3.4 nm for a = 0.35 um).

This script runs probe9 unchanged except that the lateral end-segment length is forced to --a-end (um), so that
a long device can be solved on the short-device end mesh, and vice versa a short device can be refined with
--level.  Usage (same arguments as probe9 plus --a-end):
    python phase3/c9/diag9_mesh.py --a-end 0.175 --deck R --L 1 --vg 1.32 --halo-n 2e18 --out <dir>
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from phase3.c9 import probe9 as P9  # noqa: E402  (imports devsim; installs probe9's lateral_map9 wrapper)


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument('--a-end', type=float, required=True, help='lateral graded end-segment length (um)')
    a, rest = ap.parse_known_args()
    a_end = a.a_end

    def lateral_map_fixed(L_um, a_um=0.35, h0_nm=0.5, hmax_nm=20.0):
        return P9._LMAP0(L_um, a_end, h0_nm, hmax_nm)

    P9.WP.lateral_map = lateral_map_fixed
    sys.argv = [sys.argv[0]] + rest
    P9.main()


if __name__ == '__main__':
    main()
