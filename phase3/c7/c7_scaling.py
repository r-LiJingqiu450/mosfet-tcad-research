"""C-7: empirical scaling kappa = dL_2D / l0 versus the inversion parameter alpha_inf (Ward's depletion-to-
inversion field ratio), pooled over every Ward-truncation (deck W, abrupt) state of C-6 and C-7.
The near-threshold band (3.0-3.7) was read off the C-6/C-7 structures and then tested out of sample on the
scaled device (N_A = 1e18, t_ox = 2 nm; Addendum A). Output: results/phase3/C7/analysis/{kappa_table.csv, C7_kappa.png}."""
from __future__ import annotations

import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import math
from phase3.c7.model_m7 import Device, VT

OUT = 'results/phase3/C7/analysis'


def main():
    rows = []
    for vg, dl in ((1.0, 369.873), (1.25, 153.687), (1.5, 92.878), (2.0, 53.385), (3.0, 29.909)):
        rows.append(dict(source='C6 W1 (L=4, lev0)', tox=20.0, NA=1e16, VG=vg, dL=dl, sample='development'))
    for d, smp in (('results/phase3/C7/m7', 'M7 test states'), ('results/phase3/C7/scaled', 'out-of-sample (Addendum A)'),
                   ('results/phase3/C7/kappa_test', 'out-of-sample (Addendum B)')):
        for p in sorted(glob.glob(os.path.join(d, 'e7_W_*lev0*.json'))):
            r = json.load(open(p)); a = r['args']
            rows.append(dict(source=os.path.basename(d), tox=a['tox'], NA=a['na'], VG=a['vg'],
                             dL=r['result']['const_mu']['dL_net_nm'], sample=smp))
    for z in rows:
        dev = Device(na=z['NA'], tox_nm=z['tox'])
        _, pinf, Wd = dev.far_field(z['VG'])
        l0 = dev.l0(Wd) * 1e7
        al = dev.ainf * math.sqrt(dev.ll) * dev.lam ** ((1 - pinf / (VT * dev.ll)) / 2)
        z.update(l0_nm=l0, Wd_nm=Wd * 1e7, alpha_inf=al, kappa=z['dL'] / l0, dL_over_Wd=z['dL'] / (Wd * 1e7),
                 psi_s_inf=pinf, ns_boltzmann_surface=1e10 * math.exp(pinf / VT), psi_bi_minus_psi_s=dev.psibi - pinf)
    os.makedirs(OUT, exist_ok=True)
    # kappa_table.csv = the pre-registered development table (without Addendum B); the B points go to a separate file
    dev_rows = [z for z in rows if z['sample'] != 'out-of-sample (Addendum B)']
    for fn, data in (('kappa_table.csv', dev_rows), ('kappa_table_with_B.csv', rows)):
        with open(os.path.join(OUT, fn), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(data)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.color': '#e6e5e0', 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.edgecolor': '#6b6a65'})
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    style = {(20.0, 1e16): ('#2a78d6', 'o', 't_ox 20 nm, N_A 1e16'), (10.0, 1e16): ('#eb6834', 's', 't_ox 10 nm, N_A 1e16'),
             (40.0, 1e16): ('#1baf7a', 'D', 't_ox 40 nm, N_A 1e16'), (20.0, 4e16): ('#eda100', '^', 't_ox 20 nm, N_A 4e16'),
             (2.0, 1e18): ('#1f1f1e', '*', 't_ox 2 nm, N_A 1e18 (out of sample A)'),
             (5.0, 3e17): ('#6250d6', 'P', 't_ox 5 nm, N_A 3e17 (out of sample B)'),
             (10.0, 1e17): ('#e87ba4', 'X', 't_ox 10 nm, N_A 1e17 (out of sample B)')}
    for key, (col, mk, lab) in style.items():
        e = [z for z in rows if (z['tox'], z['NA']) == key]
        ax.semilogx([z['alpha_inf'] for z in e], [z['kappa'] for z in e], mk, color=col, ms=8 if mk == '*' else 7, label=lab, ls='none')
    ax.axhspan(3.0, 3.7, xmin=0.0, xmax=1.0, color='#e6e5e0', alpha=0.5, zorder=0)
    ax.set_xlabel('α∞ (depletion-to-inversion field ratio, Ward)'); ax.set_ylabel('κ = ΔL / l₀')
    ax.set_title('Near-threshold end correction scales with l₀', loc='left', fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'C7_kappa.png'), dpi=170)
    for z in rows:
        print(f"{z['sample']:28s} tox {z['tox']:4g} NA {z['NA']:.0e} VG {z['VG']:.4f} alpha {z['alpha_inf']:6.2f} l0 {z['l0_nm']:6.1f} dL {z['dL']:6.1f} kappa {z['kappa']:.2f} dL/Wd {z['dL_over_Wd']:.2f}")


if __name__ == '__main__':
    main()
