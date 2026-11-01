"""C-12 lightweight check K8: what the simplified V_t roll-off formulas of Liang's dissertation leave out.

Setting: symmetric undoped DG, midgap gates, t_si = 10 nm, t_i = 1.5 nm SiO2, as in [L06] Figs. 3.3, 3.4 and 3.14.
All quantities are evaluated inside the Liang-Taur model [LT04 eqs. (12), (13), (16) = L06 (2.25)-(2.30)].

Constant-current roll-off as defined in [L06] eqs. (3.6) and (3.10):
    dVt(L) = -V_T ln[ (I L)(L) / (I L)(L_ref) ],
with I from LT04 eq. (16) (the full double integral).  Compared with
    (3.7)  minimum-potential method:  dVt_min = -2 sqrt(b1 c1) exp(-pi L / 2 lam1)            (symmetric, x_c = 0)
    (3.15) Taylor-expansion method:  dVt_T = dVt_min + V_T ln[ sqrt(pi)/(2 D1 L) (erf(D1 (L - y_c)) + erf(D1 y_c)) ]
and with the end-correction term alone (deduction from eq. (16), C-12 Notes N6):
    dVt_end = -V_T ln[ L / (L - dL_s - dL_d) ],   dL_{s,d} = (lam1/pi) F(b1/V_T or c1/V_T, s)   (C-12 K3)
N6: int_0^L N_inf/N dy = L - dL_s - dL_d + O(L b1 c1 / V_T^2 exp(-pi L / lam1)); the remainder is second order in the
central barrier lowering 2 sqrt(b1 c1) exp(-pi L/2 lam1) / V_T.  The first-order (minimum-potential) lowering is
already contained in the overlapping end tails, so dVt_min + dVt_end DOUBLE-COUNTS; that sum is printed only to show it.
dVt_end decays algebraically (~ 1/L); (3.7) and (3.15) decay exponentially or faster.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import erf

from phase3.c12.lt_endcorr import EPS_OX, VT, F_slice, b1_closed_form, lam1, one_minus_inv_M, m_scaled

TSI, TI, EI = 10.0, 1.5, EPS_OX


def inv_N_integral(L: float, b1: float, c1: float, l1: float, s: float) -> float:
    """int_0^L N_inf/N(y) dy for the lowest-mode LT potential (13); N(y)/N_inf = M(u_c(y)/V_T), u_c = centre value."""
    k = math.pi / l1

    def uc(y):
        # [b1 sinh(k(L-y)) + c1 sinh(k y)] / sinh(k L), overflow-safe
        a = math.exp(-k * y) * (1 - math.exp(-2 * k * (L - y)))
        b = math.exp(-k * (L - y)) * (1 - math.exp(-2 * k * y))
        return (b1 * a + c1 * b) / (1 - math.exp(-2 * k * L))

    f = lambda y: math.exp(-uc(y) / VT) / m_scaled(uc(y) / VT, s)
    pts = [p for p in (L / 2,) if 0 < p < L]
    return quad(f, 0.0, L, points=pts, limit=800, epsabs=1e-14, epsrel=1e-11)[0]


def rows_over_slice(L: float, b1: float, c1: float, l1: float, s: float) -> float:
    """I_rows / I_slice for the frozen LT potential (13) at any V_ds (Slotboom form; both share the factor
    1 - exp(-V_ds/V_T)).  Rayleigh: I_rows <= I_exact(psi) <= I_slice = LT04 eq. (16)."""
    k = math.pi / l1
    th = s * np.polynomial.legendre.leggauss(80)[0]
    w = 0.5 * np.polynomial.legendre.leggauss(80)[1]

    def uc(y):
        a = math.exp(-k * y) * (1 - math.exp(-2 * k * (L - y)))
        b = math.exp(-k * (L - y)) * (1 - math.exp(-2 * k * y))
        return (b1 * a + c1 * b) / (1 - math.exp(-2 * k * L))

    inv_rows = []
    for c in np.cos(th):
        R = quad(lambda y: math.exp(-c * uc(y) / VT), 0.0, L, points=[L / 2], limit=800, epsabs=1e-14, epsrel=1e-11)[0]
        inv_rows.append(1.0 / R)
    I_rows = float(np.sum(w * np.array(inv_rows)))          # <1/R_x>_x  (film-averaged row conductance)
    I_slice = 1.0 / inv_N_integral(L, b1, c1, l1, s)
    return I_rows / I_slice


def run(A0: float = 0.38, Lref: float = 200.0) -> dict:
    l1 = lam1(TSI, TI, EI)
    s = math.pi * TSI / (2 * l1)
    K = b1_closed_form(TSI, TI, EI, 1.0)
    out = {'tsi_nm': TSI, 'ti_nm': TI, 'eps_i': EI, 'lambda1_nm': l1, 'A0_V': A0, 'Lref_nm': Lref, 'cases': {}}
    for vds in (0.05, 1.0):
        b1, c1 = K * A0, K * (A0 + vds)
        dLs = l1 / math.pi * F_slice(b1 / VT, s)
        dLd = l1 / math.pi * F_slice(c1 / VT, s)
        IL = lambda L: L / inv_N_integral(L, b1, c1, l1, s)          # (I L) up to an L-independent factor
        ILref = IL(Lref)
        rows = []
        for L in (20.0, 24.0, 30.0, 40.0, 50.0, 70.0, 100.0, 140.0, 200.0, 400.0, 1000.0):
            full = -VT * math.log(IL(L) / ILref)
            dmin = -2 * math.sqrt(b1 * c1) * math.exp(-math.pi * L / (2 * l1))
            dmin_ref = -2 * math.sqrt(b1 * c1) * math.exp(-math.pi * Lref / (2 * l1))
            D1 = math.sqrt(math.pi ** 2 * math.sqrt(b1 * c1) * math.exp(-math.pi * L / (2 * l1)) / (l1 ** 2 * VT))
            yc = L / 2 - l1 / (2 * math.pi) * math.log(c1 / b1)
            taylor = dmin + VT * math.log(math.sqrt(math.pi) / (2 * D1 * L) * (erf(D1 * (L - yc)) + erf(D1 * yc)))
            D1r = math.sqrt(math.pi ** 2 * math.sqrt(b1 * c1) * math.exp(-math.pi * Lref / (2 * l1)) / (l1 ** 2 * VT))
            ycr = Lref / 2 - l1 / (2 * math.pi) * math.log(c1 / b1)
            taylor_ref = dmin_ref + VT * math.log(math.sqrt(math.pi) / (2 * D1r * Lref) * (erf(D1r * (Lref - ycr)) + erf(D1r * ycr)))
            endc = (-VT * math.log(L / (L - dLs - dLd)) + VT * math.log(Lref / (Lref - dLs - dLd))) if L > dLs + dLd else float('nan')
            rows.append({'L_nm': L, 'L_over_lam1': L / l1, 'dVt_full_mV': 1e3 * full,
                         'dVt_minpot_mV': 1e3 * (dmin - dmin_ref), 'dVt_taylor_mV': 1e3 * (taylor - taylor_ref),
                         'dVt_minpot_plus_endcorr_mV': 1e3 * (dmin - dmin_ref + endc),
                         'endcorr_only_mV': 1e3 * endc,
                         # exact long-channel reference: (I L)(L -> inf) = 1 in these units
                         'dVt_full_Lref_inf_mV': 1e3 * (-VT * math.log(IL(L))),
                         'I_rows_over_I_slice': rows_over_slice(L, b1, c1, l1, s)})
        out['cases'][f'Vds={vds}'] = {'b1_V': b1, 'c1_V': c1, 'dL_s_nm': dLs, 'dL_d_nm': dLd, 'rows': rows}
    return out


def main(argv=None):
    res = run()
    outdir = Path(__file__).resolve().parents[2] / 'results' / 'phase3' / 'C12'
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'lt_vt_rolloff.json').write_text(json.dumps(res, indent=1))
    print(f"lam1 = {res['lambda1_nm']:.3f} nm, A0 = {res['A0_V']} V, L_ref = {res['Lref_nm']} nm")
    for k, c in res['cases'].items():
        print(f"{k}: b1={c['b1_V']:.4f} c1={c['c1_V']:.4f} dL_s={c['dL_s_nm']:.2f} dL_d={c['dL_d_nm']:.2f} nm")
        print('   L   L/lam  full  minpot  taylor  min+end(double-counts)  endcorr | full(Lref=inf) [mV] | I_rows/I_slice')
        for r in c['rows']:
            print(f"  {r['L_nm']:6.0f} {r['L_over_lam1']:5.2f} {r['dVt_full_mV']:8.2f} {r['dVt_minpot_mV']:7.2f} {r['dVt_taylor_mV']:7.2f} "
                  f"{r['dVt_minpot_plus_endcorr_mV']:9.2f} {r['endcorr_only_mV']:8.2f} | {r['dVt_full_Lref_inf_mV']:8.2f} | "
                  f"{r['I_rows_over_I_slice']:.4f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
