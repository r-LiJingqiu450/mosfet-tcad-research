"""Figures for the C-6 W1/R1 probe. Reads results/phase3/C6/analysis_w1/w1_summary.json and the raw
w1_/r1_ JSON files. Palette: dataviz reference categorical slots 1-4 (validated; contrast relief by
distinct marker shapes, legends and the report tables)."""
from __future__ import annotations

import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from phase3.c5 import ward64_quadrature as W

C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
INK, MUTED, GRID = '#1f1f1e', '#6b6a65', '#e6e5e0'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': .6,
                     'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 2})


def main(d='results/phase3/C6/w1', a='results/phase3/C6/analysis_w1'):
    S = json.load(open(os.path.join(a, 'w1_summary.json')))
    rows = [z for z in S['rows'] if z['hw'] == 0.5]
    best = {}
    for z in rows:                       # per V_G: largest L, then finest level, with reference_ok
        if not z['reference_ok']:
            continue
        k = z['VG']; c = best.get(k)
        if c is None or (z['L'], z['level']) > (c['L'], c['level']):
            best[k] = z
    vg = np.array(sorted(best)); B = [best[v] for v in vg]
    vgf = np.linspace(1.0, 3.0, 41)
    ward = [W.dL_end_nm(W.ws_from_vg(v)) for v in vgf]

    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
    # (a) dL(V_G)
    x = ax[0]
    x.plot(vgf, ward, color=C3, label="Ward 1988 §6.4, a priori")
    x.plot(vg, [z['dL_slice'] for z in B], 's--', color=C2, ms=6, lw=1.2, label='slice (exact 2-D charge) — upper bound')
    x.plot(vg, [z['dL_Ward_hybrid'] for z in B], 'D:', color=C4, ms=6, lw=1.2, label="Ward charge + 2-D ψ_s")
    x.plot(vg, [z['dL_net'] for z in B], 'o-', color=C1, ms=7, label='2-D drift–diffusion (exact)')
    r1 = [json.load(open(p)) for p in glob.glob(os.path.join(d, 'r1_*lw0.3.json'))]
    if r1:
        x.plot([r['args']['vg'] for r in r1], [r['decomposition']['dL_net_nm'] for r in r1], 'x', color=INK, ms=8,
               mew=1.5, label='resolved n⁺ wells (R-II)', lw=0)
    x.set_xlabel('V_G (V)'); x.set_ylabel('ΔL per end (nm)'); x.set_ylim(0, None)
    x.set_title('(a) Zero-bias end correction', loc='left', color=INK)
    x.legend(frameon=False, fontsize=7.5)
    # (b) surface potential near the source
    x = ax[1]
    for v, col, mk in ((1.0, C1, 'o'), (1.25, C2, 's'), (2.0, C3, 'D')):
        p = os.path.join(d, f'w1_L4_lev0_VG{v:g}_hw0.5.json')
        if not os.path.exists(p):
            continue
        r = json.load(open(p)); s = r['decomposition']['surface']
        xs = np.array(s['x_cm']) * 1e7; ps = np.array(s['psi_V'])
        pc = ps[np.argmin(np.abs(xs - xs[-1] / 2))]
        m = xs <= 600
        x.semilogy(xs[m], (ps[m] - pc) * 1e3, '-', color=col, lw=2, label=f'2-D, V_G = {v:g} V')
        pw = (r['psi_bi_V'] - pc) * np.exp(-math.pi * xs[m] / (2 * 20.0))
        x.semilogy(xs[m], pw * 1e3, '--', color=col, lw=1.2, label=f"Ward profile, {v:g} V")
    x.set_ylim(1e-2, 3e2); x.set_xlim(0, 600)
    x.set_xlabel('distance from the n⁺ edge (nm)'); x.set_ylabel('ψ_s − ψ_s,∞ (mV)')
    x.set_title('(b) Surface potential near the source', loc='left', color=INK)
    x.legend(frameon=False, fontsize=7, ncol=1)
    # (c) GCA validity length for 1 % conductance error
    x = ax[2]
    eta = 0.01
    x.semilogy(vg, [2 * z['dL_net'] / eta / 1e3 for z in B], 'o-', color=C1, ms=7, label='2-D: L* = 2ΔL/η')
    x.semilogy(vgf, [2 * w / eta / 1e3 for w in ward], color=C3, label='Ward §6.4')
    epsL = W.EPSL_UM
    x.axhline(epsL / 0.01, color=MUTED, lw=1.2, ls='--')
    x.text(2.95, epsL / 0.01 * 1.08, 'Ward ε ≤ 0.01 criterion', ha='right', va='bottom', color=MUTED, fontsize=7.5)
    x.set_xlabel('V_G (V)'); x.set_ylabel('L* for 1 % GCA error (µm)')
    x.set_title('(c) Linear-regime GCA validity length', loc='left', color=INK)
    x.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    out = os.path.join(a, 'C6_W1_figure.png')
    fig.savefig(out, dpi=170)
    print(out)


if __name__ == '__main__':
    main(*sys.argv[1:])
