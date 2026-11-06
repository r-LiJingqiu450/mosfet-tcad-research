#!/usr/bin/env python3
"""POST-HOC EXPLORATORY DIAGNOSTIC (not pre-registered; used in no registered decision).

Question: the two G_DD forms of P1 (edge energy phi^T L phi vs -(I_S + I_B)) differ by ~1e-10 relative in the
production runs, while the two G_SD forms (-I_S vs source-excited reciprocity) agree to ~1e-14. Which form carries
the round-off?

Method (body B0, so that no current leaves through the body and every vertical cut carries the same current):
solve the drain-excited zero-bias network exactly as P1 does (benchmark/wp2/rdeck_common.terminal_solve) and evaluate
the current through every vertical cut between adjacent node lines, sum_e c_e (phi_right - phi_left) over the edges
crossing the cut. In the channel (0 <= x <= L) these differences are O(1)/N and carry no cancellation, so the cuts are
independent evaluations of G0 on the same solution. Also records the edge-conductance dynamic range.

Usage (repository root):  python benchmark/wp2/forms_diagnostic.py --id AR-1 --na 1e18 --tox 2 --L 1 --vg 0.488174 --level 0 --out DIR
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rdeck_common as RC  # noqa: E402
from rdeck_common import D, WP  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    for k in ('id', 'out'):
        ap.add_argument(f'--{k}', required=True)
    for k in ('na', 'tox', 'L', 'vg'):
        ap.add_argument(f'--{k}', type=float, required=True)
    ap.add_argument('--level', type=int, required=True)
    a = ap.parse_args()
    t0 = time.time()
    _, geo, sol0, rec = RC.build(a.L, a.level, a.tox, a.na, a.vg, 'B0')
    xs, ys = map(np.array, D.coords(D.SI))
    n0i, n1i, c, asym, psi, n = WP.network(geo)
    M = RC.terminal_masks(geo, xs, ys)
    T = dict(source=np.where(M['source'])[0], drain=np.where(M['drain'])[0])
    Lap = RC.laplacian(len(xs), n0i, n1i, c)
    rd = RC.terminal_solve(Lap, n0i, n1i, c, T, 'drain')
    phi = rd['phi']
    E, G_src = rd['energy'], -rd['I']['source']
    xl = np.asarray(geo.xs, float)
    L = geo.L_cm
    cuts = []
    for k in range(len(xl) - 1):
        m = (xs[n0i] == xl[k]) & (xs[n1i] == xl[k + 1])
        m2 = (xs[n1i] == xl[k]) & (xs[n0i] == xl[k + 1])
        I = float(np.sum(c[m] * (phi[n1i[m]] - phi[n0i[m]])) + np.sum(c[m2] * (phi[n0i[m2]] - phi[n1i[m2]])))
        xm = 0.5 * (xl[k] + xl[k + 1])
        region = 'channel' if 0.0 <= xm <= L else ('source_side' if xm < 0 else 'drain_side')
        cuts.append(dict(x_mid_um=xm * 1e4, region=region, I=I, rel_vs_source=I / G_src - 1.0, rel_vs_energy=I / E - 1.0))
    ch = [x for x in cuts if x['region'] == 'channel']
    half = [x for x in ch if x['x_mid_um'] <= 0.5 * L * 1e4]          # phi <~ 1/2: differences carry no cancellation
    out = dict(
        note='post-hoc exploratory diagnostic; not pre-registered and not used in any registered decision',
        args=vars(a), mesh=dict(nx=rec['mesh']['nx'], ny=rec['mesh']['ny']), psi_hat=rec['onecol']['psi_hat'],
        G_energy=E, G_source=G_src, energy_vs_source_rel=E / G_src - 1.0, network_relative_residual=rd['residual'],
        edge_conductance=dict(max=float(c.max()), min_nonzero=float(c[c > 0].min()), max_over_G=float(c.max() / G_src)),
        channel_cuts=dict(n=len(ch), max_abs_rel_vs_source=max(abs(x['rel_vs_source']) for x in ch),
                          median_rel_vs_source=float(np.median([x['rel_vs_source'] for x in ch])),
                          median_rel_vs_energy=float(np.median([x['rel_vs_energy'] for x in ch])),
                          max_abs_rel_vs_energy=max(abs(x['rel_vs_energy']) for x in ch)),
        source_half_channel_cuts=dict(n=len(half), max_abs_rel_vs_source=max(abs(x['rel_vs_source']) for x in half),
                                      min_rel_vs_energy=min(x['rel_vs_energy'] for x in half),
                                      max_rel_vs_energy=max(x['rel_vs_energy'] for x in half)),
        source_side_cuts_max_abs_rel_vs_source=max(abs(x['rel_vs_source']) for x in cuts if x['region'] == 'source_side'),
        drain_side_cuts_max_abs_rel_vs_source=max(abs(x['rel_vs_source']) for x in cuts if x['region'] == 'drain_side'),
        drain_side_cuts=[x for x in cuts if x['region'] == 'drain_side'][:12],
        wall_s=time.time() - t0)
    os.makedirs(a.out, exist_ok=True)
    p = os.path.join(a.out, f'forms_diag_{a.id}_B0_lev{a.level}.json')
    with open(p, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1)
    h = out['source_half_channel_cuts']
    print(f"{a.id} lev{a.level}: energy vs source {out['energy_vs_source_rel']:+.3e}; source-half channel cuts ({h['n']}) vs source max "
          f"{h['max_abs_rel_vs_source']:.1e}, vs energy [{h['min_rel_vs_energy']:+.2e}, {h['max_rel_vs_energy']:+.2e}]; channel median vs energy "
          f"{out['channel_cuts']['median_rel_vs_energy']:+.2e}; "
          f"drain-side cuts max {out['drain_side_cuts_max_abs_rel_vs_source']:.1e}; c_max/G {out['edge_conductance']['max_over_G']:.1e}", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
