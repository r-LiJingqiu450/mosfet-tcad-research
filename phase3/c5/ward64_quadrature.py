"""A-priori evaluation of Ward (1988, PhD thesis, Caltech) Sec. 6.4, eqs. (6.33)-(6.34), at v_ds -> 0.

This evaluates the PUBLISHED source-end mobile-charge formula (thesis printed p. 139) for the
F0 parameter set. It is not a derivation and not a validation; it produces the a-priori predictions
that a 2-D test (E-W1, C5 decision report) would confirm or refute.

From (6.33), ln(lam) phi_n' Q(y) = -I_ds/I_c, so at v_ds -> 0 the channel resistance is
R ~ int_0^1 dy / |Q_u(y)|. With the uniform charge (6.34) and symmetric ends at v_ds = 0,
    1 - G_GCA/G_W = 2 dL_W / L,   dL_W = (eps L) * int_0^inf (1 - |Q_m| / |Q_cs(yt)|) dyt,
where eps L = L_d sqrt(ln lam / lam) is L-independent, so dL_W is an L-independent length.

Published ingredients (p. 139):
  Q_cs(yt) ~ -sqrt(2) exp(v_s(yt) ln lam / 2) (sqrt(1 + a^2) - a),
  a(yt) = a0(yt) (ln lam)^(1/2) lam^((1 - v_s(yt))/2),
  a0(s) = sqrt(2) d_inf / (pi sqrt(s)) K((s-1)/s),  yt(s) = 2 d_inf / (pi sqrt(s)) K(1/s),  s in (1, inf),
  d_inf = sqrt(2) (2 + ln(ln lam)/ln lam)^(1/2),
  v_s(yt) = w_sl + (w_bi - w_sl) exp(-pi yt / (2 d))  [p. 139]  or  exp(-pi eps yt / (2 d))  [p. 121].
  Q_m = Q_c(w_sl : 0) = limit yt -> inf.

Interpretive choices (reported as variants; neither is stated in the thesis):
  K convention: parameter m (scipy ellipk(m)). Base choice. It reproduces Ward's stated limit
    du/dy -> -d_inf of (6.36). The modulus convention K(k) = ellipk(k^2) is reported as variant.
  Decay: the p.121 form (physical decay length 2 t_ox / pi) is the base choice. The p.139 form as
    printed (decay length 2 eps t_ox / pi, L-dependent) is reported as a variant.
  Gate bias -> w_sl: exact 1-D charge (Kingston), flat-band at V_G = 0 (F0 decks: psi_gate = V_G + psi_b0).
Asymptotic accuracy: Ward's construction carries O(1/ln lam) ~ 0.07 and O(ln ln lam / ln lam) ~ 0.19
relative errors and omits the O(1/ln lam) inner corner (6.28). Predictions are therefore order-of-
magnitude-plus (about +-30 %), not precise.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ellipk, ellipkm1

Q_E = 1.602176634e-19
VT = 0.025852
EPS_SI = 11.7 * 8.8541878128e-14
EPS_OX = 3.9 * 8.8541878128e-14
NA, NI, ND = 1e16, 1e10, 1e20
TOX_CM = 20e-7

LAM = NA / NI
LL = np.log(LAM)
LD_CM = np.sqrt(EPS_SI * VT / (Q_E * NI))
EPSL_UM = LD_CM * np.sqrt(LL / LAM) * 1e4          # eps * L in micrometres
TOX_UM = TOX_CM * 1e4
WBI = np.log(ND / NI) / LL
DINF = np.sqrt(2.0) * np.sqrt(2.0 + np.log(LL) / LL)
A_INF = DINF / np.sqrt(2.0)                          # a0 as yt -> inf (s -> 1)


def ws_from_vg(vg: float) -> float:
    """Scaled surface potential w_s (relative to intrinsic, units ln lam * V_T) from exact 1-D charge."""
    cox = EPS_OX / TOX_CM
    p0, n0 = NA, NI * NI / NA
    phif = VT * np.log(NA / NI)

    def qs(b):
        f = p0 * (np.exp(-b / VT) + b / VT - 1) + n0 * (np.exp(b / VT) - b / VT - 1)
        return np.sqrt(2 * Q_E * EPS_SI * VT) * np.sqrt(f)

    b = brentq(lambda b: b + qs(b) / cox - vg, 1e-9, 2.0)
    return (b - phif) / VT / LL


def make_conv(conv: str):
    if conv == 'parameter':
        K = ellipk
    elif conv == 'modulus':
        K = lambda k: ellipk(k * k)
    else:
        raise ValueError(conv)

    def yt_of_s(s):
        return 2 * DINF / (np.pi * np.sqrt(s)) * K(1.0 / s)

    def a0_of_s(s):
        # K((s-1)/s); for the parameter convention evaluated as ellipkm1(1/s) = K(1 - 1/s), which stays
        # finite and accurate as s -> inf (yt -> 0), where a0 -> 0 (numerical-robustness fix, C-6)
        if conv == 'parameter':
            return np.sqrt(2) * DINF / (np.pi * np.sqrt(s)) * ellipkm1(1.0 / s)
        return np.sqrt(2) * DINF / (np.pi * np.sqrt(s)) * K((s - 1.0) / s)

    def a0_of_yt(yt):
        ls = brentq(lambda ls: yt_of_s(1 + np.exp(ls)) - yt, -60, 80)
        return a0_of_s(1 + np.exp(ls))
    return a0_of_yt


def qabs(v, a0):
    a = a0 * np.sqrt(LL) * LAM ** ((1 - v) / 2)
    return np.sqrt(2) * np.exp(v * LL / 2) * (np.sqrt(1 + a * a) - a)


def dL_end_nm(ws: float, conv='parameter', decay='p121', L_um=None, mech='both') -> float:
    a0_of_yt = make_conv(conv)
    if decay == 'p121':
        dec = 2 * TOX_UM / np.pi / EPSL_UM                     # in yt units
    elif decay == 'p139_literal':
        eps = EPSL_UM / L_um
        dec = 2 * eps * TOX_UM / np.pi / EPSL_UM
    else:
        raise ValueError(decay)
    qm = qabs(ws, A_INF)

    def f(yt):
        v = ws + (WBI - ws) * np.exp(-yt / dec) if mech in ('both', 'v') else ws
        a0 = a0_of_yt(yt) if mech in ('both', 'alpha') else A_INF
        return 1 - qm / qabs(v, a0)
    val = quad(f, 1e-10, 60, limit=500, points=[dec, 1.0, DINF])[0]
    return EPSL_UM * val * 1e3


def main(out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    rec = dict(constants=dict(lam=LAM, ln_lam=LL, L_d_um=LD_CM * 1e4, epsL_um=EPSL_UM, w_bi=WBI,
                              d_inf=DINF, t_ox_um=TOX_UM, N_D=ND),
               note='A-priori evaluation of Ward 1988 Sec. 6.4 (6.33)-(6.34) at v_ds -> 0; not a validation.',
               points=[])
    for vg in (1.0, 1.25, 1.5, 2.0):
        ws = ws_from_vg(vg)
        p = dict(V_G=vg, ws_lnlam=ws * LL, inversion=bool(ws > 1.0))
        if ws > 1.0:
            p['dL_end_nm'] = dict(
                base=dL_end_nm(ws),
                modulus_K=dL_end_nm(ws, conv='modulus'),
                p139_literal_L1um=dL_end_nm(ws, decay='p139_literal', L_um=1.0),
                p139_literal_L10um=dL_end_nm(ws, decay='p139_literal', L_um=10.0),
                fringe_only=dL_end_nm(ws, mech='v'),
                charge_sharing_only=dL_end_nm(ws, mech='alpha'))
            b = p['dL_end_nm']['base']
            p['predicted_1_minus_GGCA_over_G2D'] = {f'L={L}um': 2 * b / (L * 1e3) for L in (1, 2, 4, 8, 10)}
        rec['points'].append(p)
    # check of Ward (6.36) small-x asymptote quoted on p.143 (affected by missing p.142)
    from scipy.special import ellipe
    s = 1e4
    x = (2 / np.pi) * (2 * np.sqrt(s) * (ellipe(1 / s) - ellipk(1 / s)) + ellipk(1 / s) / np.sqrt(s))
    rec['check_6_36_large_s'] = dict(s=s, x_over_dinf_times_s1p5=x * s ** 1.5,
                                     ward_quoted=-3 / 8, parameter_convention_value=-1 / 8)
    with open(os.path.join(out_dir, 'ward64_predictions.json'), 'w') as fh:
        json.dump(rec, fh, indent=2)
    print(json.dumps(rec, indent=2))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'results/phase3/C5/ward64')
