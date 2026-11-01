"""C-11 evaluation of the pre-registered predictors H-lin and H-q2D (phase3/specs/C11_Preregistration_2026-10-31.md)
on the 24 deep-subthreshold states of C-10, plus the descriptive far-tail amplitude B_inf of each 2-D solution.
Writes results/phase3/C11/c11_eval.json.
"""
from __future__ import annotations

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np  # noqa: E402
from scipy.integrate import quad  # noqa: E402

from phase3.c11 import bpredict as BP  # noqa: E402
from phase3.c7.model_m7 import Device, VT, ESI  # noqa: E402

OUT = 'results/phase3/C11'


def q2d_profile(dev, vg, top, n=4000):
    """Model M7 step 3 (C-7) with an explicit junction-face value 'top' (relative to intrinsic)."""
    psig, pinf, Wd = dev.far_field(vg)
    l0 = dev.l0(Wd)
    Cd = ESI / Wd
    Weff = l0 ** 2 * (dev.cox + Cd) / ESI
    R = lambda u: dev.Qs(u) - dev.cox * (psig - u)
    d = np.geomspace(top - pinf, 1e-12, n)
    ps = pinf + d
    Fs = np.zeros(n); acc = 0.0
    for i in range(n - 1, -1, -1):
        lo = pinf if i == n - 1 else ps[i + 1]
        acc += quad(R, lo, ps[i], limit=50)[0]
        Fs[i] = acc
    dxdpsi = np.sqrt(ESI * Weff / (2 * np.maximum(Fs, 1e-300)))
    x = np.zeros(n); x[1:] = np.cumsum(0.5 * (dxdpsi[1:] + dxdpsi[:-1]) * (ps[:-1] - ps[1:]))
    return x, ps - pinf, pinf, Wd, l0


def far_tail_B(state, lam0_nm):
    """Asymptotic amplitude from the far tail 1e-6 < u < 1e-4 V; decay length fitted there (reported)."""
    r = json.load(open(state['file']))
    x = np.array(r['result']['surface']['x_cm']) * 1e7; p = np.array(r['result']['surface']['psi_V'])
    m = (x >= 0) & (x <= state['L'] * 500 + 1e-9); x, p = x[m], p[m]
    u = p - p[-1]
    sel = (u > 1e-6) & (u < 1e-4)
    if sel.sum() < 3:
        return None
    A = np.polyfit(x[sel], np.log(u[sel]), 1)
    lam_far = -1 / A[0]
    return dict(lam_far_nm=lam_far, B_inf_V=math.exp(A[1]), B_inf_at_lam0_V=float(np.median(u[sel] * np.exp(x[sel] / lam0_nm))))


def main():
    os.makedirs(OUT, exist_ok=True)
    R = json.load(open('results/phase3/C10/analysis/c10_result.json'))
    test = [s for s in R['all_states'] if s['alpha'] >= 100]
    assert len(test) == 24, len(test)
    rows = []
    for s in test:
        dev = Device(s['na'], s['tox'])
        psig, pinf, Wd = dev.far_field(s['vg'])
        b_s = (pinf - dev.psib) / VT
        lin = BP.predict(s['na'], s['tox'], b_s, s['psibi_hat_eff'], s['hw_real_um'] * 1e3)
        x, u, pinf2, Wd2, l0 = q2d_profile(dev, s['vg'], dev.psibi + s['shift'])
        xs = np.linspace(0, x[-1], 40001)
        uu = np.interp(xs, x, u)
        r_ = s['r']
        dl_q2d = float(np.trapezoid(1 - np.exp(-uu / VT) * (1 - r_ * uu / VT), xs)) * 1e7
        ft = far_tail_B(s, lin['lam_nm'])
        rows.append(dict(file=s['file'], na=s['na'], tox=s['tox'], hw=s['hw'], shift=s['shift'], psi_hat=s['psi_hat'], alpha=s['alpha'],
                         dL_2D=s['dL'], lam_lin=lin['lam_nm'], B_lin_VT=lin['B_over_VT'], dL_lin=lin['dL_pred_nm'], err_lin=lin['dL_pred_nm'] / s['dL'] - 1,
                         dL_q2d=dl_q2d, err_q2d=dl_q2d / s['dL'] - 1, lam_fit=s['lam_fit_nm'], B_fit_VT=s['B_over_VT'],
                         far_tail=ft, B_inf_over_B_lin=(ft['B_inf_at_lam0_V'] / VT / lin['B_over_VT']) if ft else None))
    res = dict(rows=rows)
    for key in ('lin', 'q2d'):
        e = np.array([r[f'err_{key}'] for r in rows])
        n5 = int(np.sum(np.abs(e) <= 0.05)); n10 = int(np.sum(np.abs(e) > 0.10)); med = float(np.median(np.abs(e)))
        verdict = 'SUPPORTED' if (n5 >= 22 and n10 == 0) else 'REJECTED' if med > 0.05 else 'INCONCLUSIVE'
        res[f'H_{key}'] = dict(n_within_5pct=n5, n_outside_10pct=n10, median_abs_err=med, min_err=float(e.min()), max_err=float(e.max()), verdict=verdict)
    json.dump(res, open(os.path.join(OUT, 'c11_eval.json'), 'w'), indent=1)
    print('state | dL_2D | H-lin dL err | H-q2D dL err | lam_lin lam_fit lam_far | B_lin B_fit B_inf (V_T)')
    for r in sorted(rows, key=lambda z: (z['na'], z['tox'], z['hw'], z['shift'], z['psi_hat'])):
        ft = r['far_tail'] or {}
        print(f"{r['na']:.0e} tox{r['tox']:g} hw{r['hw']:g} sh{r['shift']:+.2f} psi{r['psi_hat']:5.1f} | {r['dL_2D']:7.2f} | {r['dL_lin']:7.2f} {100 * r['err_lin']:6.1f}% | {r['dL_q2d']:7.2f} {100 * r['err_q2d']:6.1f}% |"
              f" {r['lam_lin']:6.2f} {r['lam_fit']:6.2f} {ft.get('lam_far_nm', float('nan')):6.2f} | {r['B_lin_VT']:6.1f} {r['B_fit_VT']:6.1f} {ft.get('B_inf_at_lam0_V', float('nan')) / VT:6.1f}")
    for key in ('lin', 'q2d'):
        print(f'H-{key}', res[f'H_{key}'])


if __name__ == '__main__':
    main()
