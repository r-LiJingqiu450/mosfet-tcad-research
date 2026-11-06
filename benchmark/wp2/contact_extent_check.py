#!/usr/bin/env python3
"""Bound the effect of the source/drain contact extent on the zero-bias conductances (WP2, before production runs).

P1/P2 inherit the existing deck-R convention: the ohmic far-edge contacts end on the depth line nearest h_w = 0.5 um,
which is 0.498259-0.503363 um depending on device and mesh level. The device definition (and P4) uses exactly 0.5 um.
This check runs P4 (independent solver) on NON-production states with the contact ending at each realised P1 extent
and at 0.5 um, on one mesh that contains all of these depths as lines, and reports the relative changes.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import fv_rdeck as F  # noqa: E402

STATES = (  # (label, L_um, tox_nm, N_A, V_G, realised P1/P2 extents in um)
    ('dev1 L1 VG0.55 (non-production)', 1.0, 2.0, 1e18, 0.55, (0.498259, 0.503363)),
    ('dev2 L2 VG0.55 (non-production)', 2.0, 20.0, 1e16, 0.55, (0.501479, 0.503278)),
)


def main():
    out = []
    for label, L, tox, na, vg, ext in STATES:
        hard = tuple(e * 1e-4 for e in ext)
        ref = F.run(L, tox, na, vg, y_hard=hard)
        for e in ext:
            r = F.run(L, tox, na, vg, contact_extent_um=e, y_hard=hard)
            row = dict(state=label, contact_extent_um=e, reference_extent_um=0.5)
            for b in ('B0', 'B1'):
                for q in ('G_DD', 'G_SD') + (('G_DB',) if b == 'B1' else ()):
                    row[f'{b}_{q}_rel_change'] = r['results'][b][q] / ref['results'][b][q] - 1
            out.append(row)
            print(json.dumps(row), flush=True)
    worst = max(abs(v) for row in out for k, v in row.items() if k.endswith('rel_change') and 'G_DB' not in k)
    worst_db = max(abs(v) for row in out for k, v in row.items() if k.endswith('G_DB_rel_change'))
    res = dict(rows=out, max_abs_rel_change_G_DD_G_SD=worst, max_abs_rel_change_G_DB=worst_db)
    os.makedirs(os.path.join(ROOT, 'benchmark/wp2/results'), exist_ok=True)
    with open(os.path.join(ROOT, 'benchmark/wp2/results/contact_extent_check.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print(f'max |rel change| G_DD/G_SD = {worst:.2e}; G_DB = {worst_db:.2e}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
