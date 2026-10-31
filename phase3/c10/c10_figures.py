"""C-10 figure (results/phase3/C10/analysis/c10_figure.png/.pdf) from c10_result.json.

Panels:
(a) bias range: dL/L_D vs normalised band bending for device 1 (1e16 cm^-3, 20 nm) and device 2 (1e18, 2 nm),
    with the Fjeldly-Shur eq. (9) depletion width x_s/L_D for each;
(b) similarity at psi_hat = 21.02: dL/L_D of device 1 and of the device-2 variants (A..D) and the post-hoc n_i image;
(c) oxide scaling at fixed W_d: dL/W_d vs t_hat at psi_hat = 30.50 and 21.02, against the C-8/C-9 rule 1.5-1.6 W_d;
(d) far-field identity: kappa_fit = dL/lam_fit vs ln(B/V_T) for all states with alpha >= 100, line y = x + gamma;
(e) numerical issue E10-1: drain-side (pre-registered) vs source-side end correction;
(f) connection coefficient beta = B/(psi_bi - psi_s,inf) vs psi_hat.
"""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = 'results/phase3/C10/analysis'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
GAMMA = 0.5772156649015329


def style(ax, title, xl, yl):
    ax.set_title(title, loc='left', fontsize=9.5, color=INK)
    ax.set_xlabel(xl, fontsize=8.5, color=INK2); ax.set_ylabel(yl, fontsize=8.5, color=INK2)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2); ax.spines[s].set_linewidth(0.7)
    ax.tick_params(colors=INK2, labelsize=7.5)


def main():
    R = json.load(open(os.path.join(OUT, 'c10_result.json')))
    S = R['all_states']
    fig, axs = plt.subplots(2, 3, figsize=(11.6, 6.9), constrained_layout=True)
    fig.patch.set_facecolor('white')
    # (a)
    ax = axs[0, 0]
    for na, col, lab in (('1e+18', BLUE, 'device 2 (1e18, 2 nm)'), ('1e+16', ORANGE, 'device 1 (1e16, 20 nm)')):
        p = R['P10_4'][na]; LD = 4.0885 if na == '1e+18' else 40.8845
        ph = np.array(p['psi_hat']); ax.plot(ph, np.array(p['dL']) / LD, color=col, lw=1.3, marker='o', ms=4, label=f'{lab}: 2-D')
        ax.plot(ph, np.array(p['x_s_FS93']) / LD, color=col, lw=0.8, marker='x', ms=4, alpha=0.9, label=f'{lab}: FS93 eq. (9)')
    ax.legend(fontsize=6.6, frameon=False, loc='center right')
    style(ax, '(a) bias range and the FS93 depletion width', 'ψ̂ = (ψ_s∞ − ψ_b0)/V_T', 'ΔL / L_D')
    # (b)
    ax = axs[0, 1]
    d1 = R['P10_2']['device1_hw05']['dLLD']
    labels = ['device 1 (ĥ 12.3)', 'device 1 pair (ĥ 11.5)'] + [x['tag'].replace(' (both)', '') for x in R['P10_2']['rows']] + ['n_i image (all matched)']
    vals = [d1, R['P10_1']['rows'][0]['dLLD_1']] + [x['dLLD'] for x in R['P10_2']['rows']] + [R['P10_1_posthoc_ni_image']['rows'][0]['dLLD_image']]
    cols = [ORANGE, ORANGE, BLUE, BLUE, BLUE, BLUE, AQUA]
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=cols, height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v + 0.05, yi, f'{v:.2f} ({100 * (v / d1 - 1):+.1f} %)', va='center', fontsize=7, color=INK)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7); ax.set_xlim(9.0, 12.6)
    style(ax, '(b) similarity at ψ̂ = 21.02 (V_G = 0.7835 V)', 'ΔL / L_D', '')
    # (c)
    ax = axs[0, 2]
    for psi, col in ((30.50, BLUE), (21.02, AQUA)):
        pts = sorted([(s['t_hat'], s['dL'] / s['Wd_nm']) for s in S if s['na'] == 1e18 and s['hw'] == 0.5 and s['shift'] == 0 and s['level'] == 0
                      and abs(s['psi_hat'] - psi) < 0.02], key=lambda z: z[0])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, lw=1.3, marker='o', ms=5, label=f'2-D, ψ̂ = {psi:.2f}')
    ax.axhspan(1.43, 1.61, color=GRID, alpha=0.9, lw=0, label='C-8/C-9 rule ΔL ≈ 1.43–1.61 W_d')
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    style(ax, '(c) oxide scaling at fixed W_d', 't̂ = (ε_Si/ε_ox) t_ox / L_D   (t_ox = 1, 2, 6 nm)', 'ΔL / W_d')
    # (d)
    ax = axs[1, 0]
    fam = {}
    for s in S:
        if s['alpha'] < 100 or 'ln_B_VT' not in s:
            continue
        key = 'device 1' if s['na'] == 1e16 else ('t_ox 1, 6 nm' if s['tox'] != 2 else 'device 2 variants')
        fam.setdefault(key, []).append((s['ln_B_VT'], s['kappa_fit']))
    for key, col in (('device 1', ORANGE), ('device 2 variants', BLUE), ('t_ox 1, 6 nm', AQUA)):
        if key in fam:
            a = np.array(fam[key]); ax.plot(a[:, 0], a[:, 1], ls='none', marker='o', ms=5, color=col, label=key)
    xx = np.linspace(2.7, 4.6, 10)
    ax.plot(xx, xx + GAMMA, color=INK2, lw=1.0, label='κ = ln(B/V_T) + γ')
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    style(ax, '(d) far-field identity (α∞ ≥ 100)', 'ln(B / V_T)  (tail fit)', 'κ_fit = ΔL / λ_fit')
    # (e)
    ax = axs[1, 1]
    pts = sorted([(s['alpha'], s['dL_drain_prereg'] / s['dL'], s['dL_rows'] / s['dL'], s['dL_slice'] / s['dL']) for s in S if s['level'] == 0
                  and s['na'] == 1e18 and s['hw'] == 0.5 and s['shift'] == 0 and s['tox'] == 2])
    a = np.array(pts)
    ax.fill_between(a[:, 0], a[:, 2], a[:, 3], color=GRID, alpha=0.9, lw=0, label='Rayleigh bracket [rows, slice] / ΔL')
    ax.plot(a[:, 0], a[:, 1], color=ORANGE, lw=1.2, marker='o', ms=4, label='drain-side current (pre-registered ΔL)')
    ax.plot(a[:, 0], np.ones(len(a)), color=BLUE, lw=1.2, marker='o', ms=3, label='source-side current (corrected)')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.legend(fontsize=6.8, frameon=False, loc='upper left')
    style(ax, '(e) E10-1: cancellation in the drain-side current', 'α∞ (depletion/inversion field ratio)', 'ΔL / ΔL_source-side')
    # (f)
    ax = axs[1, 2]
    for key, cond, col in (('device 1', lambda s: s['na'] == 1e16, ORANGE), ('device 2 (ĥ 123)', lambda s: s['na'] == 1e18 and s['tox'] == 2 and s['hw'] == 0.5 and s['shift'] == 0, BLUE),
                           ('device 2, other ĥ/ψ̂_bi', lambda s: s['na'] == 1e18 and s['tox'] == 2 and (s['hw'] != 0.5 or s['shift'] != 0), YELLOW),
                           ('t_ox 1, 6 nm', lambda s: s['na'] == 1e18 and s['tox'] != 2, AQUA)):
        pts = sorted([(s['psi_hat'], s['beta']) for s in S if cond(s) and s['level'] == 0 and 'beta' in s and s['alpha'] >= 15])
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], ls='none' if 'other' in key or 't_ox' in key else '-', lw=1.0, marker='o', ms=4, color=col, label=key)
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    style(ax, '(f) connection coefficient is not universal', 'ψ̂', 'β = B / (ψ_bi − ψ_s∞)')
    fig.suptitle('C-10: zero-bias subthreshold end correction — direct 2-D vs Fjeldly–Shur (1993), similarity, oxide scaling, far-field identity',
                 fontsize=10, color=INK)
    fig.savefig(os.path.join(OUT, 'c10_figure.png'), dpi=300, facecolor='white')
    fig.savefig(os.path.join(OUT, 'c10_figure.pdf'), facecolor='white')
    print('wrote', os.path.join(OUT, 'c10_figure.png'))


if __name__ == '__main__':
    main()
