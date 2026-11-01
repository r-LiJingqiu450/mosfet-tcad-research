"""C-11 figure: results/phase3/C11/c11_figure.png/.pdf
(a) local decay length of the 2-D surface-potential tail vs u/V_T (device 2, psi_hat 12 and 30.5) against the 1-D
    eigen-length lam_lin: the conductance-relevant zone u ~ V_T is pre-asymptotic;
(b) parameter-free predictions (H-lin, H-q2D) vs the 2-D end correction for the 24 test states, +-5 % band;
(c) independent finite-volume solver vs the DEVSIM-based values (net, slice, rows) at psi_hat 21 and 12.
"""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = 'results/phase3/C11'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
VT = 0.025851999786435316


def style(ax, title, xl, yl):
    ax.set_title(title, loc='left', fontsize=9.5, color=INK)
    ax.set_xlabel(xl, fontsize=8.5, color=INK2); ax.set_ylabel(yl, fontsize=8.5, color=INK2)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2); ax.spines[s].set_linewidth(0.7)
    ax.tick_params(colors=INK2, labelsize=7.5)


def local_lambda(path):
    r = json.load(open(path)); x = np.array(r['result']['surface']['x_cm']) * 1e7; p = np.array(r['result']['surface']['psi_V'])
    m = x <= r['args']['L'] * 500 + 1e-9; x, p = x[m], p[m]; u = p - p[-1]
    out = []
    for i in range(len(x) - 2):
        if u[i] > 1e-6 and u[i + 2] > 1e-7 and u[i] < 0.6:
            out.append((u[i] / VT, -(x[i + 2] - x[i]) / np.log(u[i + 2] / u[i])))
    return np.array(out)


def main():
    E = json.load(open(os.path.join(OUT, 'c11_eval.json')))
    fig, axs = plt.subplots(1, 3, figsize=(11.6, 3.7), constrained_layout=True)
    fig.patch.set_facecolor('white')
    ax = axs[0]
    for f, col, lab, lam in (('results/phase3/C10/runs/e7_W_L1_lev0_VG0.488174_tox2_NA1e+18_hw0.5_ov0_sg0_dpsi0.json', BLUE, 'ψ̂ = 12', None),
                             ('results/phase3/C8/subth/e7_W_L1_lev0_VG1.08_tox2_NA1e+18_hw0.5_ov0_sg0_dpsi0.json', ORANGE, 'ψ̂ = 30.5', None)):
        a = local_lambda(f)
        ax.plot(a[:, 0], a[:, 1], color=col, lw=1.3, marker='o', ms=3, label=f'2-D tail, {lab}')
    for r in E['rows']:
        if r['na'] == 1e18 and r['tox'] == 2 and r['hw'] == 0.5 and r['shift'] == 0 and abs(r['psi_hat'] - 12) < 0.1:
            ax.axhline(r['lam_lin'], color=BLUE, lw=0.8)
        if r['na'] == 1e18 and r['tox'] == 2 and r['hw'] == 0.5 and r['shift'] == 0 and abs(r['psi_hat'] - 30.5) < 0.1:
            ax.axhline(r['lam_lin'], color=ORANGE, lw=0.8)
    ax.axvspan(0.3, 3, color=GRID, alpha=0.6, lw=0, label='zone that sets ΔL (u ≈ V_T)')
    ax.set_xscale('log'); ax.set_ylim(5, 35)
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    style(ax, '(a) the tail is pre-asymptotic where it matters', 'u / V_T  (surface-potential excess)', 'local decay length (nm)   [lines: λ_lin]')
    ax = axs[1]
    d = np.array([r['dL_2D'] for r in E['rows']])
    for key, col, lab in (('dL_lin', BLUE, 'H-lin (eigen-projection)'), ('dL_q2d', ORANGE, 'H-q2D (quasi-2-D, M7)')):
        p = np.array([r[key] for r in E['rows']])
        ax.plot(d, p / d, ls='none', marker='o', ms=5, color=col, label=f"{lab}: {E['H_' + key[3:]]['verdict']}")
    ax.axhspan(0.95, 1.05, color=GRID, alpha=0.9, lw=0, label='±5 % (registered)')
    ax.set_xscale('log'); ax.legend(fontsize=7, frameon=False, loc='lower right'); ax.set_ylim(0.6, 1.1)
    style(ax, '(b) parameter-free predictions (24 states)', 'ΔL from 2-D (nm)', 'ΔL_pred / ΔL_2D')
    ax = axs[2]
    I = {vg: json.load(open(f'{OUT}/indep/fd2d_VG{vg}.json')) for vg in ('0.783500', '0.488174')}
    ref = {'0.783500': dict(net=47.324, slice=47.327, rows=47.190), '0.488174': dict(net=41.700, slice=43.010, rows=40.749)}
    cats = ['net', 'slice', 'rows']; xx = np.arange(3)
    for k, (vg, col, lab) in enumerate((('0.783500', AQUA, 'ψ̂ = 21'), ('0.488174', BLUE, 'ψ̂ = 12'))):
        rel = [100 * (I[vg][f'dL_{c}_nm'] / ref[vg][c] - 1) for c in cats]
        ax.bar(xx + (k - 0.5) * 0.36, rel, width=0.34, color=col, label=lab)
    ax.axhline(0, color=INK2, lw=0.8); ax.set_xticks(xx); ax.set_xticklabels(['network G', 'slice bound', 'rows bound'], fontsize=8)
    ax.legend(fontsize=7, frameon=False)
    style(ax, '(c) independent solver vs DEVSIM-based ΔL', '', 'difference (%)')
    fig.suptitle('C-11: the subthreshold end correction is not predicted by linear or quasi-2-D reduced models', fontsize=10, color=INK)
    fig.savefig(os.path.join(OUT, 'c11_figure.png'), dpi=300, facecolor='white'); fig.savefig(os.path.join(OUT, 'c11_figure.pdf'), facecolor='white')
    print('wrote', os.path.join(OUT, 'c11_figure.png'))


if __name__ == '__main__':
    main()
