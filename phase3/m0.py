"""M0: semi-analytic GCA (Pao-Sah first-integral) reference with injected parameters.

Package B reconstruction (2026-10-30), Research_05 v2.1 B-2. Rebuilt from the
Phase 2C source because the 2026-10-29 Package B files were not preserved.

This module re-implements, without modification of semantics,
``semiconductor_sim/stage3_mosfet/reference.py`` @ 7073f20:
  * ``_em1mx``            -> reference.py lines 16-19
  * ``reference_slice``   -> reference.py lines 22-44
  * ``reference_current`` -> reference.py lines 47-56
Physical constants and defaults come from ``common/constants.py`` (read-only).
The only change is that every physical parameter and every tolerance is passed
explicitly (``M0Params``, ``Tolerances``) instead of read from module globals.
With ``M0Params.f0()`` and ``Tolerances()`` the arithmetic is the same sequence
of operations as production, so results should agree bitwise (test B-2).

Units: SI internally (m, m^-3, V, A). ``I`` is the current for width W (A);
``I'`` = I / W[um] is reported in A/um. Uncertainty numbers are estimates, not
bounds (Research_04 v2.1 section 7.1).
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, replace

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

import phase3  # noqa: F401  (sets sys.path for production imports)
from common import constants as C


@dataclass(frozen=True)
class M0Params:
    na: float          # acceptor density (m^-3)
    tox: float         # oxide thickness (m)
    x_max: float       # body depth (m)
    L: float           # channel length (m)
    W: float           # channel width (m)
    mu_n: float        # electron mobility (m^2/Vs)
    vfb: float         # flat-band voltage (V)
    ni: float = C.NI
    vt: float = C.VT
    eps_si: float = C.EPS_SI
    eps_ox: float = C.EPS_OX
    q: float = C.Q_E

    @property
    def cox(self) -> float:
        return self.eps_ox / self.tox

    @property
    def p0(self) -> float:
        return .5 * (self.na + np.sqrt(self.na * self.na + 4 * self.ni ** 2))

    @property
    def psi_b0(self) -> float:
        """Bulk potential referenced to the intrinsic level (V); F0: -0.357159 V."""
        return -self.vt * math.log(self.p0 / self.ni)

    @classmethod
    def f0(cls) -> "M0Params":
        """Family F0 = production defaults (constants.py, reference.py x_max=1.2 um)."""
        return cls(na=C.NA, tox=C.EPS_OX / C.COX, x_max=1.2e-6, L=C.L_CH,
                   W=C.W_CH, mu_n=C.MU_N, vfb=C.VFB)

    def to_dict(self):
        d = asdict(self)
        d.update(cox=self.cox, p0=self.p0, psi_b0=self.psi_b0)
        return d


@dataclass(frozen=True)
class Tolerances:
    """Production tolerances (reference.py). ``tightened(f)`` divides all by f."""
    brentq_xtol: float = 2e-14
    slice_epsabs: float = 1e-6
    slice_epsrel: float = 2e-10
    slice_limit: int = 200
    cur_excess_epsabs: float = 1e-4
    cur_bg_epsabs: float = 1e-6
    cur_epsrel: float = 2e-8
    cur_limit: int = 160

    def tightened(self, f: float = 100.) -> "Tolerances":
        return replace(self, slice_epsabs=self.slice_epsabs / f,
                       slice_epsrel=self.slice_epsrel / f,
                       cur_excess_epsabs=self.cur_excess_epsabs / f,
                       cur_bg_epsabs=self.cur_bg_epsabs / f,
                       cur_epsrel=self.cur_epsrel / f,
                       slice_limit=self.slice_limit * 4, cur_limit=self.cur_limit * 4)


DEFAULT_TOL = Tolerances()


def _em1mx(z):
    """exp(z)-1-z with a Taylor branch (reference.py lines 16-19)."""
    if abs(z) < 1e-3:
        return z * z * (.5 + z * (1 / 6 + z * (1 / 24 + z * (1 / 120 + z / 720))))
    return np.expm1(z) - z


def reference_slice(vg, phi_n, p: M0Params, tol: Tolerances = DEFAULT_TOL):
    """One GCA slice by the Poisson first integral (reference.py lines 22-44).

    Potential phi is referenced to the neutral bulk (phi_bulk ~ 0).
    Returns phi_s, qs (C/m^2), ns (m^-2, incl. background nb*x_max), ns_excess,
    ns_background, nb, quadrature_abs_estimate, bulk_shift.
    """
    if vg < .05 or phi_n < 0:
        raise ValueError('Reference supports positive inversion/depletion gate bias')
    na, ni, VT, q, eps = p.na, p.ni, p.vt, p.q, p.eps_si
    p0 = .5 * (na + np.sqrt(na * na + 4 * ni ** 2))
    product = ni ** 2 * np.exp(-phi_n / VT)
    pb = .5 * (na + np.sqrt(na * na + 4 * product)); nb = product / pb
    bulk = -VT * np.log(pb / p0)

    def energy(u):
        return VT * (pb * _em1mx(-u) + nb * _em1mx(u))

    def charge(ps):
        u = (ps - bulk) / VT
        return -np.sqrt(max(0., 2 * eps * q * energy(u)))

    ps = brentq(lambda v: p.vfb + v - charge(v) / p.cox - vg,
                bulk + 1e-14, max(2.5, phi_n + 1), xtol=tol.brentq_xtol)
    us = (ps - bulk) / VT

    def sheet_integrand(u):
        if abs(u) < 1e-8:
            return nb * np.sqrt(eps * VT / (q * (pb + nb)))
        return nb * np.expm1(u) * VT / np.sqrt(2 * q / eps * energy(u))

    excess, error = quad(sheet_integrand, 0, us, epsabs=tol.slice_epsabs,
                         epsrel=tol.slice_epsrel, limit=tol.slice_limit)
    return dict(phi_s=ps, qs=charge(ps), ns=excess + nb * p.x_max, ns_excess=excess,
                ns_background=nb * p.x_max, nb=nb,
                quadrature_abs_estimate=error, bulk_shift=bulk)


def reference_current(vg, vd, p: M0Params, tol: Tolerances = DEFAULT_TOL):
    """I = (q W mu_n / L) int_0^VD Ns(VG, phi_n) dphi_n (reference.py lines 47-56).

    Returns I in A for width W. ``quadrature_abs_estimate`` covers only the
    outer integral, as in production; inner errors are probed by B-3.
    """
    excess, error = quad(lambda pn: reference_slice(vg, pn, p, tol)['ns_excess'], 0, vd,
                         epsabs=tol.cur_excess_epsabs, epsrel=tol.cur_epsrel,
                         limit=tol.cur_limit)
    background, bg_error = quad(lambda pn: reference_slice(vg, pn, p, tol)['ns_background'],
                                0, vd, epsabs=tol.cur_bg_epsabs, epsrel=tol.cur_epsrel,
                                limit=tol.cur_limit)
    pref = p.q * p.W * p.mu_n / p.L
    return dict(ID=pref * (excess + background), ID_excess=pref * excess,
                ID_background=pref * background,
                background_fraction=background / (excess + background),
                quadrature_abs_estimate=pref * (error + bg_error))


# ---------------------------------------------------------------- helpers
def w_um(p: M0Params) -> float:
    return p.W * 1e6


def current_per_um(vg, vd, p: M0Params, tol: Tolerances = DEFAULT_TOL) -> float:
    """I' in A/um. Zero at VD = 0 exactly."""
    if vd == 0:
        return 0.
    return reference_current(vg, vd, p, tol)['ID'] / w_um(p)


def conductance0(vg, p: M0Params, tol: Tolerances = DEFAULT_TOL) -> float:
    """G'_M0 = dI'/dVD at VD=0 = q mu_n Ns(VG,0) / L, in S/um (R04 v2.1 5.2(C))."""
    ns = reference_slice(vg, 0., p, tol)['ns']
    return p.q * p.W * p.mu_n * ns / p.L / w_um(p)


def ith_per_um(p: M0Params) -> float:
    """Constant-current criterion I'_th(L) = 1e-7 A / L[um] in A/um (R04 3.2)."""
    return 1e-7 / (p.L * 1e6)


def vt_cc_direct(vd, p: M0Params, lo=.3, hi=2., tol: Tolerances = DEFAULT_TOL,
                 xtol=1e-12):
    """Bracketed root of ln I'(VG) = ln I'_th on M0 itself. Returns (VT, status)."""
    target = math.log(ith_per_um(p))
    f = lambda v: math.log(current_per_um(v, vd, p, tol)) - target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return float('nan'), 'not_bracketed'
    return brentq(f, lo, hi, xtol=xtol, rtol=4 * np.finfo(float).eps), 'ok'


def vt_cc_tabulated(vg, i_per_um, ith):
    """Common finder for tabulated curves (also for 2-D data): monotone cubic (PCHIP)
    interpolation of ln I' vs VG, then a bracketed root. Returns (VT, status)."""
    vg = np.asarray(vg, float); li = np.log(np.asarray(i_per_um, float))
    if not np.all(np.diff(vg) > 0):
        return float('nan'), 'vg_not_increasing'
    if not np.all(np.diff(li) > 0):
        return float('nan'), 'current_not_monotone'
    t = math.log(ith)
    if not (li[0] < t < li[-1]):
        return float('nan'), 'not_bracketed'
    k = int(np.searchsorted(li, t))
    f = PchipInterpolator(vg, li)
    return brentq(lambda v: float(f(v)) - t, vg[k - 1], vg[k], xtol=1e-14), 'ok'


def delta_vg(vg_a, i_a, vg_b, i_b, i_target):
    """delta V_G = V_G,A(I') - V_G,B(I') at fixed current (R04 6), in V."""
    va, sa = vt_cc_tabulated(vg_a, i_a, i_target)
    vb, sb = vt_cc_tabulated(vg_b, i_b, i_target)
    return va - vb, (sa if sa != 'ok' else sb)
