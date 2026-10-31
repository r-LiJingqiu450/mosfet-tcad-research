"""C-9 figure: the four pre-registered tests (H9-1..H9-4), the split-C-V audit and the halo realism check, on the
direct 2-D solves of results/phase3/C9/runs (level 0) plus the post-hoc level-1 mesh diagnostic of
results/phase3/C9/diag/level1 where present.

(a) H9-1  affine length law: dL per end vs V_G, deck R without halo, L = 0.5 / 1 / 2 um.
(b) H9-3  weak halo: dL vs V_G at L = 1 um, no halo vs 5e17 halo (zero line shows the sign is not reversed).
(c)       strong halo (2e18): dL vs V_G at L = 1 um, level 0 (production) and level 1 (mesh diagnostic).
(d) H9-4  maximum-g_m threshold shift vs 1/L with Lombardi mobility (no halo, 5e17 halo), and the C-8 definition
          (constant mobility, extrapolation at the sweep edge) for comparison; thin lines are c/L fits.
(e)       strong halo: dV_T vs 1/L, direct points and the affine-law identity G(L) = sigma_inf/(L - 2 dL) evaluated
          from the L = 1 um solves (level 0 and level 1).
(f) H9-2  split-C-V mobility error vs V_G at L = 1 um: terminal charge vs electron proxy (decks W and R), and the
          two-length (1 and 2 um) differential.
Output: results/phase3/C9/analysis/e9_figure.png (300 dpi) and e9_figure.pdf.  Data: e9_result_v2.json and runs_v2
(halo states re-run with probe9 v2 after erratum E9-1; no-halo states from runs/), and diag/level1_v2 for the mesh check.
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from phase3.c9 import e9_analysis as E  # noqa: E402

OUT = 'results/phase3/C9/analysis'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'   # categorical slots 1-4 (validated)
RAMP = {0.5: '#86b6ef', 1.0: '#2a78d6', 2.0: '#104281'}                    # ordinal blue ramp (validated)


HALO_V2 = 'results/phase3/C9/runs_v2'
LEVEL1 = 'results/phase3/C9/diag/level1_v2'


def load_level1():
    D = {}
    for p in glob.glob(os.path.join(LEVEL1, 'e9_*.json')):
        r = json.load(open(p))
        if 'end_utc' not in r:
            continue
        a = r['args']; c, m = r['result']['const_mu'], r['result']['lombardi_mu']
        D.setdefault((a['deck'], a['halo_n'], a['L']), {})[a['vg']] = dict(dL=c['dL_net_nm'], dLmu=m['dL_net_nm'], Gmu=m['G_net'],
                                                                        GGCAmu=m['G_GCA_1COL_mu'])
    return D


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc='left', fontsize=9.5, color=INK)
    ax.set_xlabel(xlabel, fontsize=8.5, color=INK2); ax.set_ylabel(ylabel, fontsize=8.5, color=INK2)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2); ax.spines[s].set_linewidth(0.7)
    ax.tick_params(colors=INK2, labelsize=7.5)


def vt_affine(v, GG1, dL1_nm, L_um, vmax=2.0):
    """V_T(L) from the affine identity using the L = 1 um solve: G(L) = sigma_inf / (L - 2 dL)."""
    return E.max_gm_vt(v, GG1 / (L_um - 2e-3 * dL1_nm), vmax)[0]


def main():
    D = E.load(HALO_V2); D1 = load_level1()
    res = json.load(open(os.path.join(OUT, 'e9_result_v2.json')))
    fig, axs = plt.subplots(2, 3, figsize=(11.0, 6.6), constrained_layout=True)
    fig.patch.set_facecolor('white')

    # (a) affine law
    ax = axs[0, 0]
    for L, mk, ms in ((0.5, 'o', 7), (1.0, 's', 4.5), (2.0, '.', 4)):
        v, dl = E.series(D, ('R', 0.0, L), 'dL')
        ax.plot(v, dl, color=RAMP[L], lw=1.2, marker=mk, ms=ms, mfc='white' if L == 0.5 else RAMP[L], label=f'L = {L:g} µm')
    dmax = {dk: max(abs(a - b) for a, b in zip(res['H9_1'][dk]['dL_05'], res['H9_1'][dk]['dL_2'])) for dk in ('R', 'W')}
    ax.text(0.97, 0.95, f"max |ΔL(0.5) − ΔL(2)| over all V_G:\n{dmax['R']:.2f} nm (deck R), {dmax['W']:.2f} nm (deck W)", transform=ax.transAxes,
            ha='right', va='top', fontsize=7.5, color=INK2)
    ax.legend(fontsize=7.5, frameon=False, loc='lower left')
    style(ax, '(a) H9-1: affine length law (no halo)', 'V_G (V)', 'ΔL per end (nm)')

    # (b) weak halo
    ax = axs[0, 1]
    for hn, col, lab in ((0.0, BLUE, 'no halo'), (5e17, ORANGE, 'halo 5e17 cm⁻³')):
        v, dl = E.series(D, ('R', hn, 1.0), 'dL')
        ax.plot(v, dl, color=col, lw=1.2, marker='o', ms=3.5, label=lab)
        ax.annotate(lab, (v[-1], dl[-1]), xytext=(4, 0), textcoords='offset points', fontsize=7.5, color=INK, va='center')
    ax.axhline(0, color=INK2, lw=0.8)
    v, dl = E.series(D, ('R', 5e17, 1.0), 'dL'); k = list(v).index(1.32)
    ax.annotate(f'+{dl[k]:.1f} nm at 1.32 V:\nsign not reversed\n(H9-3 rejected)', (1.32, dl[k]), xytext=(1.45, 36), fontsize=7.5,
                color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.7))
    ax.set_xlim(1.05, 2.25); ax.set_ylim(-5, 55)
    style(ax, '(b) H9-3: weak halo, L = 1 µm', 'V_G (V)', 'ΔL per end (nm)')

    # (c) strong halo, mesh levels
    ax = axs[0, 2]
    v, dl = E.series(D, ('R', 2e18, 1.0), 'dL')
    ax.plot(v, dl, color=AQUA, lw=1.2, marker='o', ms=3.5, label='halo 2e18, level 0')
    if ('R', 2e18, 1.0) in D1:
        v1 = sorted(D1[('R', 2e18, 1.0)]); y1 = [D1[('R', 2e18, 1.0)][x]['dL'] for x in v1]
        ax.plot(v1, y1, ls='none', marker='o', ms=7, mfc='none', mec=AQUA, mew=1.0, label='level 1 (mesh check)')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.legend(fontsize=7.5, frameon=False, loc='lower right')
    style(ax, '(c) strong halo, L = 1 µm (constant μ)', 'V_G (V)', 'ΔL per end (nm)')

    # (d) max-gm dV_T vs 1/L
    ax = axs[1, 0]
    xs = np.linspace(0, 2.15, 50)
    for key, col, lab in (('R_halo0e+00_vmax2.0', BLUE, 'no halo (Lombardi μ, max-g_m)'), ('R_halo5e+17_vmax2.0', ORANGE, 'halo 5e17')):
        rows = res['H9_4'][key]['rows']; x = np.array([1 / r['L'] for r in rows]); y = np.array([r['dVT_mV'] for r in rows])
        c = float(np.sum(x * y) / np.sum(x * x))
        ax.plot(xs, c * xs, color=col, lw=0.8); ax.plot(x, y, ls='none', marker='o', ms=5, color=col, label=lab)
    c8 = json.load(open('results/phase3/C8/extraction/e8_2_result.json'))['vt_linear_extrapolation']
    c8 = [r for r in c8 if r['L_um'] in (0.5, 1.0, 2.0)]
    ax.plot([1 / r['L_um'] for r in c8], [r['dVT_mV'] for r in c8], ls='none', marker='x', ms=6, color=INK2,
            label='C-8 definition (constant μ, sweep edge)')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.legend(fontsize=7, frameon=False, loc='lower left')
    ax.set_xlim(0, 2.15)
    style(ax, '(d) H9-4: threshold shift, no / weak halo', '1/L (µm⁻¹)', 'ΔV_T (mV)')

    # (e) strong halo dV_T vs 1/L: direct points and affine identity from L = 1 um
    ax = axs[1, 1]
    rows = res['H9_4']['R_halo2e+18_vmax2.0']['rows']
    v, Gm = E.series(D, ('R', 2e18, 1.0), 'GGCAmu'); _, dlm = E.series(D, ('R', 2e18, 1.0), 'dLmu')
    vinf = E.max_gm_vt(v, Gm, 2.0)[0]
    Ls = np.linspace(0.3, 20.0, 400)
    ax.plot(1 / Ls, [1e3 * (vt_affine(v, Gm, dlm, L) - vinf) for L in Ls], color=AQUA, lw=1.0, label='affine identity, L = 1 µm solves (level 0)')
    if ('R', 2e18, 1.0) in D1 and len(D1[('R', 2e18, 1.0)]) == len(v):
        v1 = np.array(sorted(D1[('R', 2e18, 1.0)])); G1 = np.array([D1[('R', 2e18, 1.0)][x]['GGCAmu'] for x in v1])
        d1 = np.array([D1[('R', 2e18, 1.0)][x]['dLmu'] for x in v1]); vinf1 = E.max_gm_vt(v1, G1, 2.0)[0]
        ax.plot(1 / Ls, [1e3 * (vt_affine(v1, G1, d1, L) - vinf1) for L in Ls], color=INK2, lw=1.0, label='same, level-1 solves')
    ax.plot([1 / r['L'] for r in rows], [r['dVT_mV'] for r in rows], ls='none', marker='o', ms=5, color=AQUA, label='direct 2-D (level 0), Lombardi μ')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.legend(fontsize=7, frameon=False, loc='lower right')
    ax.set_xlim(0, 2.15)
    style(ax, '(e) strong halo 2e18: reverse short-channel', '1/L (µm⁻¹)', 'ΔV_T (mV)')

    # (f) split-C-V error
    ax = axs[1, 2]
    for deck, col in (('W', YELLOW), ('R', BLUE)):
        e = res['H9_2'][f'{deck}_L1']
        ax.plot(e['VG'], 100 * np.array(e['err_terminal']), color=col, lw=1.2, marker='o', ms=4, label=f'deck {deck}: terminal charge')
        ax.plot(e['VG'], 100 * np.array(e['err_electrons']), ls='none', marker='o', ms=7, mfc='none', mec=col, mew=1.0,
                label=f'deck {deck}: electron proxy (C-8)')
    t = res['two_length_R']
    ax.plot(t['VG'], 100 * np.array(t['rel_err']), color=INK2, lw=1.0, marker='^', ms=4, label='two-length differential (R, 1 & 2 µm)')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.legend(fontsize=7, frameon=False, loc='lower right')
    style(ax, '(f) H9-2: split-C–V mobility error, L = 1 µm', 'V_G (V)', 'μ error (%)')

    fig.suptitle('C-9: direct 2-D audit of the C-8 extraction claims (N_A = 1e18 cm⁻³, t_ox = 2 nm, zero bias, F0 physics)',
                 fontsize=10, color=INK)
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, 'e9_figure.png'), dpi=300, facecolor='white')
    fig.savefig(os.path.join(OUT, 'e9_figure.pdf'), facecolor='white')
    print('wrote', os.path.join(OUT, 'e9_figure.png'))


if __name__ == '__main__':
    main()
