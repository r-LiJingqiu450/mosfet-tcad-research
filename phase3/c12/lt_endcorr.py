"""C-12 lightweight checks: the zero-drain-bias subthreshold end correction of an undoped symmetric double-gate (DG)
MOSFET inside the Liang-Taur (LT) analytical model.

Sources (full text read):
  [LT04] X. Liang and Y. Taur, "A 2-D analytical solution for SCEs in DG MOSFETs," IEEE Trans. Electron Devices,
         vol. 51, no. 8, pp. 1385-1391, Aug. 2004.  Eqs. (1)-(16) and g1(x) at the bottom of p. 1388.
  [L06]  X. Liang, "Analytical modeling of short channel effects in double gate MOSFET," Ph.D. dissertation,
         Univ. California, San Diego, 2006.  Eqs. (2.1)-(2.30); Sec. 4.1-4.4.

Everything below stays inside the LT boundary-value problem (Laplace's equation in Si and insulator, abrupt Dirichlet
source/drain faces at Eg/2q and Eg/2q + V_ds across the film, linear interpolation in the insulator gaps, gates at
V_g - dphi).  No device simulation is run.  Lengths in nm, potentials in V, symmetric DG, V_ds -> 0.

Checks
  K1  The LT "conjugate function" g1 (p. 1388) equals eps(x) u_L1(x,0) / eps_si, i.e. the eps-weighted eigenfunction;
      LT04's "non-orthogonal eigenfunctions + conjugate set" (eqs. 9-11) and L06's "orthogonal with weight eps(x)"
      (eq. 2.23) are the same projection.  Also reproduces the lambda_1 values quoted in L06 Figs. 4.2(c) and 4.5.
  K2  The closed form b1, LT04 eq. (12) = L06 eq. (2.25), equals the projection eq. (11) of the LT boundary data.
  K3  Long-channel end correction implied by eqs. (13) + (16) at V_ds -> 0 (slice form):
          dL = (lam1/pi) F(beta, s),  F = int_0^beta [1 - 1/M(z)] dz/z,  M(z) = (1/2s) int_{-s}^{s} exp(z cos th) dth,
          beta = b1/V_T,  s = pi t_si / (2 lam1).
      Thin-film limit s -> 0: F -> Ein(beta), which is the C-11 pre-registered DG predictor lam_DG Ein(B_proj/V_T).
      Large beta: F = ln(beta) + C(s) + o(1), C(0) = Euler's gamma, C(s) < gamma for s > 0.
  K4  Finite-L quadrature of eq. (16) with eq. (13): converges to K3 as L/lam1 grows; for short L the conductance
      "end correction" absorbs the short-channel barrier lowering.
  K5  Truncation: full symmetric-mode series of the same boundary-value problem vs the lowest mode kept by LT.
  K6  Transport: eq. (16) is the slice (Dirichlet) upper bound of the frozen-potential conductance; the rows (Thomson)
      lower bound is (lam1/pi) <Ein(beta cos(pi x/lam1))>_x.  Both coincide in the thin-film limit.
  K7  Sensitivity to the source-edge boundary amplitude (the dominant LT error per L06 Sec. 4.2): d(dL)/d ln b1.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import exp1

VT = 0.025851999786435316          # kT/q at 300 K (same constant as C-10/C-11)
EPS_SI = 11.7
EPS_OX = 3.9
EG2 = 0.56                          # Eg/2q used by LT for the n+ source (L06 Fig. 4.2: 0.56 V)
GAMMA = 0.5772156649015329
_GL_X, _GL_W = np.polynomial.legendre.leggauss(400)


def ein(z: float) -> float:
    """Entire exponential integral Ein(z) = int_0^z (1 - e^-t)/t dt = gamma + ln z + E1(z)."""
    if z < 1e-6:
        return z - z * z / 4.0
    return GAMMA + math.log(z) + float(exp1(z))


# ---------------------------------------------------------------- eigenproblem (symmetric modes)
def _f_sym(k: float, tsi: float, ti: float, ei: float, esi: float = EPS_SI) -> float:
    """Pole-free form of LT04 eq. (8), odd n:  eps_si tan(k ti) = eps_i cot(k tsi/2),  k = pi/lambda."""
    return esi * math.sin(k * ti) * math.sin(0.5 * k * tsi) - ei * math.cos(k * ti) * math.cos(0.5 * k * tsi)


def sym_wavenumbers(tsi: float, ti: float, ei: float, nmodes: int, esi: float = EPS_SI) -> np.ndarray:
    """First nmodes roots k_n of _f_sym (k = pi/lambda_n), by sign scan + Brent."""
    dk = min(math.pi / (0.5 * tsi), math.pi / ti) / 200.0
    roots, k0, f0 = [], 1e-9, _f_sym(1e-9, tsi, ti, ei, esi)
    while len(roots) < nmodes:
        k1 = k0 + dk
        f1 = _f_sym(k1, tsi, ti, ei, esi)
        if f0 == 0.0:
            roots.append(k0)
        elif f0 * f1 < 0.0:
            roots.append(brentq(_f_sym, k0, k1, args=(tsi, ti, ei, esi), xtol=1e-15, rtol=1e-15))
        k0, f0 = k1, f1
    return np.array(roots)


def lam1(tsi: float, ti: float, ei: float, esi: float = EPS_SI) -> float:
    return math.pi / sym_wavenumbers(tsi, ti, ei, 1, esi)[0]


def _ins_amp(k: float, tsi: float, ti: float, ei: float, esi: float = EPS_SI) -> float:
    """Amplitude C of the insulator part C sin(k xi) (xi = distance from the gate electrode) of the mode cos(k x)."""
    s = math.sin(k * ti)
    if abs(s) > 1e-6:
        return math.cos(0.5 * k * tsi) / s
    return esi * math.sin(0.5 * k * tsi) / (ei * math.cos(k * ti))   # flux-continuity form


def project_lt_data(k: float, tsi: float, ti: float, ei: float, A0: float = 1.0, esi: float = EPS_SI) -> float:
    """eps-weighted projection of the LT source-face data (A0 in Si, linear 0 -> A0 from gate to interface in the
    insulator; LT04 eq. (2) and the interpolation below it; L06 (2.4), (2.7)-(2.8)) onto the mode with wavenumber k.
    Analytic integrals."""
    C = _ins_amp(k, tsi, ti, ei, esi)
    num = esi * A0 * (2.0 / k) * math.sin(0.5 * k * tsi) \
        + 2.0 * ei * A0 * C / ti * (math.sin(k * ti) / k ** 2 - ti * math.cos(k * ti) / k)
    den = esi * (0.5 * tsi + math.sin(k * tsi) / (2.0 * k)) + 2.0 * ei * C ** 2 * (0.5 * ti - math.sin(2.0 * k * ti) / (4.0 * k))
    return num / den


def b1_closed_form(tsi: float, ti: float, ei: float, A0: float = 1.0, esi: float = EPS_SI) -> float:
    """LT04 eq. (12) (= L06 eq. (2.25)) with V_ds = 0, symmetric gates; A0 = Eg/2q + dphi - V_g."""
    l1 = lam1(tsi, ti, ei, esi)
    num = 2.0 * l1 ** 2 * math.tan(math.pi * ti / l1) * math.sin(math.pi * tsi / (2.0 * l1))
    den = math.pi ** 2 * ti * (0.5 * tsi + math.sin(math.pi * tsi / l1) / math.sin(2.0 * math.pi * ti / l1) * ti)
    return num / den * A0


def b1_quadrature_eq11(tsi: float, ti: float, ei: float, A0: float = 1.0, esi: float = EPS_SI) -> float:
    """LT04 eq. (11) evaluated by Gauss-Legendre quadrature with g1(x) exactly as printed at the bottom of p. 1388 and
    u_L1(x,0) from eq. (9).  Independent route to K2."""
    l1 = lam1(tsi, ti, ei, esi)
    a = math.pi / l1

    def g1(x):
        if x < -0.5 * tsi:
            return math.sin(a * 0.5 * tsi) / math.cos(a * ti) * math.sin(a * (x + ti + 0.5 * tsi))
        if x <= 0.5 * tsi:
            return math.cos(a * x)
        return -math.sin(a * 0.5 * tsi) / math.cos(a * ti) * math.sin(a * (x - ti - 0.5 * tsi))

    def u1(x):  # eq. (9), y = 0, n = 1
        if x < -0.5 * tsi:
            return math.sin(math.pi / 2 - a * 0.5 * tsi) / math.sin(a * ti) * math.sin(a * (x + ti + 0.5 * tsi))
        if x <= 0.5 * tsi:
            return math.sin(math.pi / 2 + a * x)
        return -math.sin(math.pi / 2 + a * 0.5 * tsi) / math.sin(a * ti) * math.sin(a * (x - ti - 0.5 * tsi))

    def data(x):  # psi(x,0) - v(x)
        if x < -0.5 * tsi:
            return A0 * (x + 0.5 * tsi + ti) / ti
        if x <= 0.5 * tsi:
            return A0
        return A0 * (0.5 * tsi + ti - x) / ti

    num = den = 0.0
    for lo, hi in ((-0.5 * tsi - ti, -0.5 * tsi), (-0.5 * tsi, 0.5 * tsi), (0.5 * tsi, 0.5 * tsi + ti)):
        xs = 0.5 * (hi - lo) * _GL_X + 0.5 * (hi + lo)
        w = 0.5 * (hi - lo) * _GL_W
        num += sum(wi * data(xi) * g1(xi) for xi, wi in zip(xs, w))
        den += sum(wi * u1(xi) * g1(xi) for xi, wi in zip(xs, w))
    return num / den


def conjugate_identity_residual(tsi: float, ti: float, ei: float, esi: float = EPS_SI) -> float:
    """max | g1(x) - eps(x) u_L1(x,0)/eps_si | / max|g1| over the insulator (K1); zero by eq. (8)."""
    l1 = lam1(tsi, ti, ei, esi)
    a = math.pi / l1
    xi = np.linspace(0.0, ti, 101)
    g = math.sin(a * 0.5 * tsi) / math.cos(a * ti) * np.sin(a * xi)
    u = math.cos(a * 0.5 * tsi) / math.sin(a * ti) * np.sin(a * xi)
    return float(np.max(np.abs(g - (ei / esi) * u)) / max(1.0, float(np.max(np.abs(g)))))


# ---------------------------------------------------------------- end-correction integrals
_TH_X, _TH_W = np.polynomial.legendre.leggauss(200)


def m_scaled(z: float, s: float) -> float:
    """m(z) = M(z) e^{-z} = (1/2s) int_{-s}^{s} exp(z (cos th - 1)) dth  (film-averaged exp of the lowest mode)."""
    if s < 1e-12:
        return 1.0
    th = s * _TH_X
    return float(0.5 * np.sum(_TH_W * np.exp(z * (np.cos(th) - 1.0))))


def one_minus_inv_M(z: float, s: float) -> float:
    if z < 1e-8:
        th = s * _TH_X
        return float(z * 0.5 * np.sum(_TH_W * np.cos(th))) if s > 1e-12 else z
    return 1.0 - math.exp(-z) / m_scaled(z, s)


def F_slice(beta: float, s: float) -> float:
    """F(beta, s) = int_0^beta [1 - 1/M(z)] dz/z  (dL_slice = (lam1/pi) F)."""
    f = lambda z: one_minus_inv_M(z, s) / z if z > 0 else (0.5 * np.sum(_TH_W * np.cos(s * _TH_X)) if s > 1e-12 else 1.0)
    pts = [p for p in (1.0, 5.0, 20.0) if p < beta]
    val, err = quad(f, 0.0, beta, points=pts or None, limit=400, epsabs=1e-13, epsrel=1e-12)
    return val


def C_const(s: float) -> float:
    """C(s) = lim_{beta->inf} [F(beta,s) - ln beta] = int_0^1 [1-1/M] dz/z - int_1^inf dz/(z M)."""
    a = quad(lambda z: one_minus_inv_M(z, s) / z if z > 0 else 1.0, 0.0, 1.0, epsabs=1e-14, epsrel=1e-12, limit=200)[0]
    b = quad(lambda z: math.exp(-z) / (m_scaled(z, s) * z), 1.0, np.inf, epsabs=1e-14, epsrel=1e-12, limit=400)[0]
    return a - b


def F_rows(beta: float, s: float) -> float:
    """Rows (Thomson) lower bound, long channel: <Ein(beta cos th)>_th over |th| < s."""
    if s < 1e-12:
        return ein(beta)
    th = s * _TH_X
    return float(0.5 * np.sum(_TH_W * np.array([ein(beta * math.cos(t)) for t in th])))


def dL_finite_L(beta: float, s: float, L_over_lam: float) -> float:
    """dL(L) = (1/2) int_0^L [1 - 1/M(beta h(y))] dy for eq. (13) at V_ds = 0, in units of lam1/pi.
    h(y) = [sinh(k(L-y)) + sinh(k y)] / sinh(k L) with k y = t (t in [0, pi L/lam1])."""
    T = math.pi * L_over_lam

    def h(t):
        # numerically stable for large T
        return (math.exp(-t) * (1 - math.exp(-2 * (T - t))) + math.exp(-(T - t)) * (1 - math.exp(-2 * t))) / (1 - math.exp(-2 * T))

    f = lambda t: one_minus_inv_M(beta * h(t), s)
    val = quad(f, 0.0, 0.5 * T, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    return val                      # = (1/2) * 2 * half-channel integral, in units of 1/k = lam1/pi


def dL_rows_finite_L(beta: float, s: float, L_over_lam: float) -> float:
    """Rows bound for finite L (units lam1/pi): dL = (L/2)(1 - G1D/Grows), G_rows/G1D = <L / int_0^L e^{-u/VT} dy>."""
    T = math.pi * L_over_lam
    th = s * _TH_X if s > 1e-12 else np.array([0.0])
    w = 0.5 * _TH_W if s > 1e-12 else np.array([1.0])

    def h(t):
        return (math.exp(-t) * (1 - math.exp(-2 * (T - t))) + math.exp(-(T - t)) * (1 - math.exp(-2 * t))) / (1 - math.exp(-2 * T))

    ratios = []
    for c in np.cos(th):
        R = 2.0 * quad(lambda t: math.exp(-beta * c * h(t)), 0.0, 0.5 * T, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
        ratios.append(T / R)
    g = float(np.sum(w * np.array(ratios)))
    return 0.5 * T * (1.0 - 1.0 / g)


def dL_full_series(tsi: float, ti: float, ei: float, A0: float, nmodes: int = 300, esi: float = EPS_SI) -> dict:
    """K5: long-channel dL (slice) with the full symmetric-mode series of the LT boundary-value problem, in nm."""
    ks = sym_wavenumbers(tsi, ti, ei, nmodes, esi)
    bs = np.array([project_lt_data(k, tsi, ti, ei, A0, esi) for k in ks])
    xs = 0.5 * tsi * _GL_X
    wx = 0.5 * _GL_W                          # weights for the film average (1/t_si) int dx
    cosk = np.cos(np.outer(ks, xs))           # (nmodes, nx)

    def one_minus(y):
        u = (bs * np.exp(-ks * y)) @ cosk / VT       # u(x, y)/V_T on the Si quadrature nodes
        um = float(np.max(u))
        avg = float(np.sum(wx * np.exp(u - um)))
        return 1.0 - math.exp(-um) / avg

    k1 = ks[0]
    ymax = (math.log(max(bs[0] / VT, 2.0)) + 40.0) / k1
    val = quad(one_minus, 0.0, ymax, limit=800, epsabs=1e-12, epsrel=1e-11, points=[1.0 / k1, 3.0 / k1])[0]
    # residual of the boundary data at y = 0 (series reproduces A0 in Si away from the corners)
    u0 = (bs @ cosk)
    return {'dL_nm': val, 'b1': float(bs[0]), 'b3_over_b1': float(bs[1] / bs[0]), 'b5_over_b1': float(bs[2] / bs[0]),
            'series_at_y0_centre_over_A0': float(u0[len(u0) // 2] / A0), 'nmodes': nmodes}


def k9_source_edge_estimate() -> dict:
    """K9 (rough estimate, not a measurement): effect of the source-edge idealisation on b1 and dL for the case of
    L06 Fig. 4.2(c), p. 78 (t_si = 5 nm, t_i = 1.5 nm SiO2, V_g = 0.2 V midgap, L = 18 nm, V_ds = 1 V).
    Values read by eye from the figure (uncertainty about +-0.01 V): LT data 0.56 V flat in Si; DESSIS 0.545 V at the
    centre and 0.49 V at the Si/insulator interfaces, roughly linear to the gate value 0.2 V across the insulator.
    The DESSIS profile in Si is modelled as a parabola.  Projection by Gauss-Legendre quadrature with the eps weight."""
    tsi, ti, ei = 5.0, 1.5, EPS_OX
    l1 = lam1(tsi, ti, ei)
    k = math.pi / l1
    C = _ins_amp(k, tsi, ti, ei)
    vg = 0.2

    def proj(f_si_centre, f_si_edge):
        num = den = 0.0
        xs = 0.5 * tsi * _GL_X
        w = 0.5 * tsi * _GL_W
        prof = (f_si_centre - (f_si_centre - f_si_edge) * (2 * xs / tsi) ** 2) - vg
        num += EPS_SI * float(np.sum(w * prof * np.cos(k * xs)))
        den += EPS_SI * float(np.sum(w * np.cos(k * xs) ** 2))
        xi = 0.5 * ti * (_GL_X + 1.0)                         # distance from the gate
        wi = 0.5 * ti * _GL_W
        prof_i = (f_si_edge - vg) * xi / ti
        num += 2 * ei * C * float(np.sum(wi * prof_i * np.sin(k * xi)))
        den += 2 * ei * C ** 2 * float(np.sum(wi * np.sin(k * xi) ** 2))
        return num / den

    b_lt = proj(EG2, EG2)
    b_sim = proj(0.545, 0.49)
    s = math.pi * tsi / (2 * l1)
    beta = b_lt / VT
    return {'lambda1_nm': l1, 'b1_LT_V': b_lt, 'b1_from_read_DESSIS_profile_V': b_sim, 'ratio': b_sim / b_lt,
            'dL_rel_change_long_channel': F_slice(beta * b_sim / b_lt, s) / F_slice(beta, s) - 1.0,
            'note': 'figure values read by eye; V_ds = 1 V and L/lam1 ~ 1.5 in the source figure; order-of-magnitude only'}


# ---------------------------------------------------------------- run
GEOMS = {
    'G1_tsi10_ti1.5_SiO2': (10.0, 1.5, EPS_OX),     # LT04 Figs. 4-7 class (t_i = 1.5 nm, t_si = 10 nm)
    'G2_tsi5_ti1.5_SiO2': (5.0, 1.5, EPS_OX),       # L06 Fig. 4.2(c): lambda1 = 11.7 nm quoted
    'G3_tsi25_ti1.5_SiO2': (25.0, 1.5, EPS_OX),     # L06 Fig. 4.5(a): lambda1 = 33.6 nm quoted
    'G4_tsi25_ti3_k7.6': (25.0, 3.0, 7.6),          # L06 Fig. 4.5(b): lambda1 = 33.9 nm quoted
    'G5_tsi25_ti13.5_k35.1': (25.0, 13.5, 35.1),    # L06 Fig. 4.5(c): lambda1 = 39.0 nm quoted
}
QUOTED_LAMBDA1 = {'G2_tsi5_ti1.5_SiO2': 11.7, 'G3_tsi25_ti1.5_SiO2': 33.6, 'G4_tsi25_ti3_k7.6': 33.9,
                  'G5_tsi25_ti13.5_k35.1': 39.0}
BARRIERS = (0.30, 0.40, 0.50)       # A0 = Eg/2q + dphi - V_g (V): deep subthreshold, as in LT's 0.2 V-below-V_t criterion


def run() -> dict:
    out = {'constants': {'VT': VT, 'eps_si': EPS_SI, 'eps_ox': EPS_OX, 'Eg_over_2q': EG2}, 'geometries': {}}
    for name, (tsi, ti, ei) in GEOMS.items():
        l1 = lam1(tsi, ti, ei)
        s = math.pi * tsi / (2.0 * l1)
        g = {'tsi_nm': tsi, 'ti_nm': ti, 'eps_i': ei, 'lambda1_nm': l1, 'lambda1_quoted_nm': QUOTED_LAMBDA1.get(name),
             'ti_over_tsi': ti / tsi, 'LT_interpolation_validity_ti_over_tsi_lt_0.3': ti / tsi < 0.3,
             's_rad': s, 'K1_conjugate_identity_residual': conjugate_identity_residual(tsi, ti, ei)}
        K_cf = b1_closed_form(tsi, ti, ei, 1.0)
        g['K2_b1_per_volt_closed_form_eq12'] = K_cf
        g['K2_b1_per_volt_projection'] = project_lt_data(math.pi / l1, tsi, ti, ei, 1.0)
        g['K2_b1_per_volt_quadrature_eq11'] = b1_quadrature_eq11(tsi, ti, ei, 1.0)
        g['K2_rel_diff_projection'] = g['K2_b1_per_volt_projection'] / K_cf - 1.0
        g['K2_rel_diff_quadrature'] = g['K2_b1_per_volt_quadrature_eq11'] / K_cf - 1.0
        g['C_s'] = C_const(s)
        g['C_s_minus_gamma'] = g['C_s'] - GAMMA
        rows = []
        for A0 in BARRIERS:
            b1 = K_cf * A0
            beta = b1 / VT
            Fs = F_slice(beta, s)
            Fr = F_rows(beta, s)
            Fe = ein(beta)
            r = {'A0_V': A0, 'b1_V': b1, 'beta': beta,
                 'dL_slice_LT_nm': l1 / math.pi * Fs, 'dL_rows_nm': l1 / math.pi * Fr, 'dL_Ein_C11_nm': l1 / math.pi * Fe,
                 'Ein_vs_LT_rel': Fe / Fs - 1.0, 'rows_slice_bracket_rel': 1.0 - Fr / Fs,
                 'asymptotic_lnbeta_plus_C_rel': (math.log(beta) + g['C_s']) / Fs - 1.0}
            # K4 finite L
            r['K4_dL_over_dLinf'] = {str(Lr): dL_finite_L(beta, s, Lr) / Fs for Lr in (3.0, 5.0, 10.0, 20.0)}
            r['K6_rows_finite_L20_over_rows_inf'] = dL_rows_finite_L(beta, s, 20.0) / Fr
            # K7 sensitivity to b1
            r['K7_rel_change_dL_for_b1_x0.8_0.9_1.1'] = [F_slice(f * beta, s) / Fs - 1.0 for f in (0.8, 0.9, 1.1)]
            r['K7_dlnDL_dlnb1'] = (F_slice(1.01 * beta, s) - F_slice(0.99 * beta, s)) / (0.02 * Fs)
            # K5 full series
            fs = dL_full_series(tsi, ti, ei, A0)
            r['K5_full_series'] = fs
            r['K5_full_over_one_mode_rel'] = fs['dL_nm'] / (l1 / math.pi * Fs) - 1.0
            rows.append(r)
        g['bias'] = rows
        out['geometries'][name] = g
    # thin-film limit of F -> Ein (K3), independent of geometry
    out['K3_thin_film_limit'] = {str(s): F_slice(15.0, s) / ein(15.0) - 1.0 for s in (1e-3, 0.05, 0.2, 0.5, 0.8, 1.2)}
    out['K3_C_of_s'] = {str(s): C_const(s) for s in (1e-3, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5)}
    out['K9_source_edge_estimate'] = k9_source_edge_estimate()
    return out


def main(argv=None):
    res = run()
    outdir = Path(__file__).resolve().parents[2] / 'results' / 'phase3' / 'C12'
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'lt_endcorr.json').write_text(json.dumps(res, indent=1))
    for name, g in res['geometries'].items():
        print(f"{name}: lam1={g['lambda1_nm']:.3f} nm (quoted {g['lambda1_quoted_nm']}), s={g['s_rad']:.3f}, "
              f"K1 res={g['K1_conjugate_identity_residual']:.1e}, b1/A0={g['K2_b1_per_volt_closed_form_eq12']:.5f} "
              f"(proj {g['K2_rel_diff_projection']:+.1e}, quad {g['K2_rel_diff_quadrature']:+.1e}), C(s)-gamma={g['C_s_minus_gamma']:+.4f}")
        for r in g['bias']:
            print(f"   A0={r['A0_V']:.2f} beta={r['beta']:.2f} dL_LT={r['dL_slice_LT_nm']:.3f} rows={r['dL_rows_nm']:.3f} "
                  f"Ein={r['dL_Ein_C11_nm']:.3f} nm | Ein-LT {100*r['Ein_vs_LT_rel']:+.2f}% bracket {100*r['rows_slice_bracket_rel']:.2f}% "
                  f"| full/1mode {100*r['K5_full_over_one_mode_rel']:+.3f}% | dlnDL/dlnb1 {r['K7_dlnDL_dlnb1']:.3f} "
                  f"| L/lam=3,5,10,20: " + ' '.join(f"{v:.4f}" for v in r['K4_dL_over_dLinf'].values()))
    print('thin-film limit F/Ein-1 at beta=15:', {k: f"{v:+.2e}" for k, v in res['K3_thin_film_limit'].items()})
    print('C(s):', {k: round(v, 5) for k, v in res['K3_C_of_s'].items()})
    k9 = res['K9_source_edge_estimate']
    print(f"K9 (estimate from L06 Fig. 4.2(c)): b1_sim/b1_LT = {k9['ratio']:.3f}, long-channel dL change {100*k9['dL_rel_change_long_channel']:+.2f}%")
    return 0


if __name__ == '__main__':
    sys.exit(main())
