"""E8-2: consequences of the intrinsic gate-dependent end correction for standard extraction methods.

Input: scaled-device sweep (N_A = 1e18, t_ox = 2 nm, deck W, abrupt, L = 1 um, V_G = 1.20..2.00 V, constant mu).
By the exact affine law (C-6 P6) the zero-bias conductance per unit width at mask/metallurgical length L is
    G(L, V_G) = sigma_inf(V_G) / (L - 2 dL(V_G)),   sigma_inf = q mu N_inf (1COL),  no series resistance.
(a) Terada-Muta: intersection of 1/G-vs-L lines at V1 < V2 -> (L*, R*); exact closed form
      R_sd,TM = 2[dL(V1) - dL(V2)] / [sigma(V2) - sigma(V1)],  dL_TM = L*/2.
(b) maximum-g_m linear-extrapolation threshold V_T(L) (cubic spline of G(V_G) on the sweep grid).
(c) split-C-V mobility: measured gate-modulated charge Q(L,V_G) = q{[N_inf(V_G)-N_inf(V_lo)] L + 2[E(V_G)-E(V_lo)]},
      E = per-end excess sheet charge int_0^{L/2} (N(x)-N_inf) dx of the 2-D solution (spill-over cancels in the
      difference); mu_split = L^2 G / Q.  True mu = 1350 cm^2/Vs everywhere.
Outputs results/phase3/C8/extraction/{e8_2_result.json, e8_2_figure.png}."""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from scipy.interpolate import CubicSpline

Q, MU = 1.602176634e-19, 1350.0
OUT = 'results/phase3/C8/extraction'


def load():
    rows = []
    for p in glob.glob('results/phase3/C8/sweep/e7_W_*.json') + glob.glob('results/phase3/C8/subth/e7_W_L1_*NA1e+18*.json'):
        r = json.load(open(p))
        c = r['result']['const_mu']; cs = r['result']['column_sheet_density']
        x = np.array(cs['x_cm']); N = np.array(cs['N_cm2']); L = x[-1]
        h = x <= L / 2 + 1e-15
        Ninf = float(N[np.argmin(np.abs(x - L / 2))])
        E = float(np.trapezoid(N[h] - Ninf, x[h]))                      # cm^-1 (per unit width), per end
        rows.append(dict(VG=r['args']['vg'], dL_cm=c['dL_net_nm'] * 1e-7, dL_nm=c['dL_net_nm'], Ns=r['result']['Ns_cm2'],
                         sigma=Q * MU * r['result']['Ns_cm2'], E_end=E, Ninf_centre=Ninf))
    rows.sort(key=lambda z: z['VG'])
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    R = load()
    vg = np.array([z['VG'] for z in R]); dL = np.array([z['dL_cm'] for z in R]); sg = np.array([z['sigma'] for z in R])
    Lset_um = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 1e6]
    res = dict(sweep=[dict(VG=z['VG'], dL_nm=z['dL_nm'], Ns=z['Ns'], E_end=z['E_end']) for z in R])
    # (a) Terada-Muta
    tm = []
    for v1, v2 in ((1.40, 1.60), (1.60, 2.00), (1.40, 2.00), (1.32, 1.52)):   # nearest sweep points
        i, j = int(np.argmin(abs(vg - v1))), int(np.argmin(abs(vg - v2)))
        Ls = 2 * (dL[i] * sg[j] - dL[j] * sg[i]) / (sg[j] - sg[i])
        Rs = 2 * (dL[i] - dL[j]) / (sg[j] - sg[i])                         # ohm*cm (per unit width)
        tm.append(dict(V1=float(vg[i]), V2=float(vg[j]), dL_true_V1_nm=dL[i] * 1e7, dL_true_V2_nm=dL[j] * 1e7,
                       dL_TM_nm=Ls / 2 * 1e7, Rsd_TM_ohm_um=Rs * 1e4, Rsq_V1_ohm=1 / sg[i], Rsq_V2_ohm=1 / sg[j],
                       Rsd_TM_over_channel_at_L1um=Rs / ((1e-4 - 2 * dL[i]) / sg[i])))
    res['terada_muta'] = tm
    # (b) maximum-gm linear-extrapolation V_T(L)
    vt = []
    for Lum in Lset_um:
        L = Lum * 1e-4
        G = sg / (L - 2 * dL) * (L if Lum > 1e5 else 1.0) if Lum > 1e5 else sg / (L - 2 * dL)
        if Lum > 1e5:
            G = sg.copy()                                                   # infinitely long: sigma only (scale-free)
        cs = CubicSpline(vg, G)
        f = np.linspace(vg[0], vg[-1], 4001); g = cs(f); gm = cs(f, 1)
        k = int(np.argmax(gm))
        vt.append(dict(L_um=Lum, VT_lin=float(f[k] - g[k] / gm[k]), VG_at_gm_max=float(f[k])))
    VTinf = vt[-1]['VT_lin']
    for z in vt:
        z['dVT_mV'] = (z['VT_lin'] - VTinf) * 1e3
        z['dVT_times_L_mV_um'] = z['dVT_mV'] * z['L_um'] if z['L_um'] < 1e5 else None
    res['vt_linear_extrapolation'] = vt
    # (c) split-C-V mobility (V_lo = lowest sweep point)
    Ninf = np.array([z['Ns'] for z in R]); E = np.array([z['E_end'] for z in R])
    sc = []
    for Lum in (0.5, 1.0, 2.0, 5.0, 10.0):
        L = Lum * 1e-4
        for k in range(1, len(R)):
            Qm = Q * ((Ninf[k] - Ninf[0]) * L + 2 * (E[k] - E[0]))
            G = sg[k] / (L - 2 * dL[k])
            mu = L * L * G / Qm
            mu_long = sg[k] / (Q * (Ninf[k] - Ninf[0]))          # same procedure, infinitely long device
            sc.append(dict(L_um=Lum, VG=float(vg[k]), mu_split=mu, rel_err=mu / MU - 1, rel_err_end_only=mu / mu_long - 1,
                           dL_C_nm=float((E[k] - E[0]) / (Ninf[k] - Ninf[0]) * 1e7), dL_nm=float(dL[k] * 1e7)))
    res['split_cv'] = sc
    json.dump(res, open(os.path.join(OUT, 'e8_2_result.json'), 'w'), indent=1)
    print('Terada-Muta:'); [print({k: round(v, 3) for k, v in t.items()}) for t in tm]
    print('V_T(L):'); [print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in z.items()}) for z in vt]
    print('split-CV (L, VG, rel_err, dL_C, dL):')
    for z in sc:
        if z['VG'] in (1.32, 1.4, 1.6, 2.0) or abs(z['VG'] - 1.32) < 1e-9:
            print(z['L_um'], z['VG'], round(z['rel_err'], 4), 'end-only', round(z['rel_err_end_only'], 4), round(z['dL_C_nm'], 1), round(z['dL_nm'], 1))
    # figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.color': '#e6e5e0', 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.5))
    ax[0].plot(vg, dL * 1e7, 'o-', color='#2a78d6', ms=4, label='ΔL (conductance)')
    dLC = [(E[k] - E[0]) / (Ninf[k] - Ninf[0]) * 1e7 for k in range(1, len(R))]
    ax[0].plot(vg[1:], dLC, 's--', color='#eb6834', ms=4, label='ΔL_C (gate-modulated charge)')
    ax[0].set_xlim(1.18, 2.02); ax[0].set_yscale('log')
    ax[0].set_xlabel('V_G (V)'); ax[0].set_ylabel('per-end length (nm)'); ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title('(a) End lengths, N_A 1e18, t_ox 2 nm', loc='left', fontsize=9)
    Ls = np.array([z['L_um'] for z in vt[:-1]]); dv = np.array([z['dVT_mV'] for z in vt[:-1]])
    ax[1].plot(1 / Ls, dv, 'o-', color='#2a78d6', ms=5)
    ax[1].set_xlabel('1 / L (µm⁻¹)'); ax[1].set_ylabel('ΔV_T,lin(L) (mV)')
    ax[1].set_title('(b) Max-g_m threshold shift from end term only', loc='left', fontsize=9)
    for Lum, col in ((0.5, '#2a78d6'), (1.0, '#eb6834'), (2.0, '#1baf7a'), (10.0, '#eda100')):
        e = [z for z in sc if z['L_um'] == Lum and z['VG'] >= 1.28]
        ax[2].plot([z['VG'] for z in e], [100 * z['rel_err_end_only'] for z in e], color=col, label=f'L = {Lum:g} µm')
    ax[2].set_xlabel('V_G (V)  [N_s: 1.4e11 at 1.32 V, 5.5e11 at 1.40 V]'); ax[2].set_ylabel('split-C–V mobility error (%)'); ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].set_title('(c) Split-C–V μ: end-term error vs long device', loc='left', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'e8_2_figure.png'), dpi=170)


if __name__ == '__main__':
    main()
