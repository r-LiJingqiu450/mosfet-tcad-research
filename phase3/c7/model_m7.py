"""Model M7 (C-7): parameter-free reduced model of the zero-bias source/drain end correction dL(V_G)
of a long-channel MOSFET with abrupt n+ junctions (Ward truncation), F0 physics.

Status: WORKING HYPOTHESIS. Every ingredient is either exact 1-D electrostatics, an exact linear
eigenvalue of the depletion-approximation end problem, or Ward's published corner-charge formula
([W88] p. 139, Sec. 6.3). No parameter is fitted to 2-D data.

  1. Far field: exact 1-D charge Q_s(psi) (Boltzmann, frozen holes), C_ox (psi_g - psi_inf) = Q_s(psi_inf),
     psi_g = V_G + psi_b0 (V_FB = 0 convention of the decks).
  2. Natural length l0 of the end layer in the depletion approximation (lowest decaying mode of the
     two-layer strip: oxide (eps_ox, t_ox, gate Dirichlet) over depleted Si (eps_si, W_d, delta-psi = 0 at the
     depletion edge to first order)):  tan(k W_d) = -(eps_si/eps_ox) tan(k t_ox),  k in (pi/2W_d, pi/W_d), l0 = 1/k.
     W_d = 1-D depletion width at psi_inf (depletion approximation incl. the -V_T term).
  3. Surface potential psi_s(x): generalised quasi-2-D equation
        eps_si W_eff psi'' = Q_s(psi) - C_ox (psi_g - psi),   psi(0) = psi_bi,  psi(inf) = psi_inf,
     with W_eff = l0^2 (C_ox + C_d)/eps_si, C_d = eps_si/W_d, so that its subthreshold linearisation decays with
     exactly l0. Solved by its first integral (psi')^2 = 2F(psi)/(eps_si W_eff).  In inversion the exact
     Q_s(psi) automatically shortens the decay (inversion screening).
  4. Local inversion charge: Ward's thin-inversion-layer corner formula with the 90-degree corner
     depletion-field factor alpha_0l(y~) ([W88] p. 139; K-parameter convention as in phase3/c5):
        |Q| ~ sqrt(2) exp(v ln(lam)/2) (sqrt(1 + a^2) - a),  a = alpha_0l (ln lam)^(1/2) lam^((1-v)/2),
     with v = psi_s(x)/(V_T ln lam) from step 3 (instead of Ward's oxide-fringe profile).
  5. dL_M7 = int_0^inf (1 - |Q_inf|/|Q(x)|) dx  (slice form, Ward (6.33)); current spreading not modelled.
Domain of intended validity (stated before testing): threshold to moderate inversion, i.e. surface
inversion charge not larger than the depletion charge (alpha_inf >~ 0.5); end layers not overlapping.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.optimize import brentq
from scipy.special import ellipk, ellipkm1

Q = 1.602176634e-19
VT = 0.025851999786435316
EPS0 = 8.8541878128e-14
ESI, EOX = 11.7 * EPS0, 3.9 * EPS0
NI = 1e10
ND = 1e20
LD_CM = math.sqrt(ESI * VT / (Q * NI))


class Device:
    def __init__(self, na=1e16, tox_nm=20.0, nd=ND):
        self.na, self.tox = na, tox_nm * 1e-7
        self.p0 = .5 * (na + math.sqrt(na * na + 4 * NI * NI)); self.n0 = NI * NI / self.p0
        self.psib = -VT * math.log(self.p0 / NI)
        self.cox = EOX / self.tox
        self.psibi = VT * math.log(nd / NI)
        self.lam = na / NI; self.ll = math.log(self.lam)
        self.epsL = LD_CM * math.sqrt(self.ll / self.lam)             # cm
        self.dinf = math.sqrt(2.0) * math.sqrt(2.0 + math.log(self.ll) / self.ll)
        self.ainf = self.dinf / math.sqrt(2.0)

    # exact 1-D semiconductor charge magnitude (C/cm^2), psi relative to intrinsic
    def Qs(self, psi):
        b = (psi - self.psib) / VT
        f = self.p0 * (math.exp(-b) + b - 1) + self.n0 * (math.expm1(b) - b)
        return math.sqrt(2 * Q * ESI * VT * max(f, 0.0))

    def far_field(self, vg):
        psig = vg + self.psib
        psi = brentq(lambda p: self.Qs(p) - self.cox * (psig - p), self.psib + 1e-9, psig)
        bb = psi - self.psib
        Wd = math.sqrt(2 * ESI * max(bb - VT, 1e-6) / (Q * self.na))
        return psig, psi, Wd

    def l0(self, Wd):
        r = ESI / EOX
        f = lambda k: math.tan(k * Wd) + r * math.tan(k * self.tox)
        k = brentq(f, (math.pi / 2) / Wd * (1 + 1e-9), math.pi / Wd * (1 - 1e-12))
        return 1.0 / k

    # Ward corner factor (parameter convention, robust as y~ -> 0)
    def a0(self, yt):
        if yt <= 0:
            return 0.0
        yts = lambda s: 2 * self.dinf / (math.pi * math.sqrt(s)) * ellipk(1.0 / s)
        if yt >= yts(1 + math.exp(-27.0)):
            return self.ainf
        ls = brentq(lambda ls: yts(1 + math.exp(ls)) - yt, -27.0, 80.0)
        s = 1 + math.exp(ls)
        return math.sqrt(2) * self.dinf / (math.pi * math.sqrt(s)) * ellipkm1(1.0 / s)

    def qward(self, psi, a0):
        v = psi / (VT * self.ll)
        a = a0 * math.sqrt(self.ll) * self.lam ** ((1 - v) / 2)
        return math.sqrt(2) * math.exp(v * self.ll / 2) * (math.sqrt(1 + a * a) - a)

    def profile(self, vg, n=4000):
        """psi_s(x) from the generalised quasi-2-D equation; returns x (cm, ascending), psi."""
        psig, pinf, Wd = self.far_field(vg)
        l0 = self.l0(Wd)
        Cd = ESI / Wd
        Weff = l0 ** 2 * (self.cox + Cd) / ESI
        R = lambda u: self.Qs(u) - self.cox * (psig - u)
        top = self.psibi
        # psi grid from psi_bi down to psi_inf + tiny, geometric in (psi - psi_inf)
        d = np.geomspace(top - pinf, 1e-9, n)
        ps = pinf + d
        # F(psi) = int_{pinf}^{psi} R(u) du, cumulative from the bottom
        Fs = np.zeros(n)
        prev = 0.0
        for i in range(n - 1, -1, -1):
            lo = pinf if i == n - 1 else ps[i + 1]
            prev += quad(R, lo, ps[i], limit=50)[0]
            Fs[i] = prev
        dxdpsi = np.sqrt(ESI * Weff / (2 * np.maximum(Fs, 1e-300)))
        x = np.concatenate([[0.0], cumulative_trapezoid(dxdpsi, -ps)]) if False else None
        x = np.zeros(n)
        x[1:] = np.cumsum(0.5 * (dxdpsi[1:] + dxdpsi[:-1]) * (ps[:-1] - ps[1:]))
        return dict(x=x, psi=ps, psi_inf=pinf, psig=psig, Wd=Wd, l0=l0, Weff=Weff)

    def dL(self, vg, xmax_factor=60.0):
        pr = self.profile(vg)
        x, ps, pinf = pr['x'], pr['psi'], pr['psi_inf']
        qinf = self.qward(pinf, self.ainf)
        # integrate on x: beyond the computed profile psi = psi_inf (deviation < 1e-9 V)
        xm = max(x[-1], xmax_factor * pr['l0'])
        xg = np.unique(np.concatenate([x, np.linspace(0, xm, 4000)]))
        pg = np.interp(xg, x, ps, right=pinf)
        a0 = np.array([self.a0(t / self.epsL) for t in xg])
        qg = np.array([self.qward(p, a) for p, a in zip(pg, a0)])
        integ = 1 - qinf / qg
        val = float(np.trapezoid(integ, xg))
        # same with the 2-D-free 'local 1-D' charge (no alpha reduction): diagnostic
        q1 = np.array([self.qward(p, self.ainf) for p in pg])
        val_psi_only = float(np.trapezoid(1 - qinf / q1, xg))
        alpha_inf = self.ainf * math.sqrt(self.ll) * self.lam ** ((1 - pinf / (VT * self.ll)) / 2)
        return dict(dL_nm=val * 1e7, dL_psi_only_nm=val_psi_only * 1e7, l0_nm=pr['l0'] * 1e7, Wd_nm=pr['Wd'] * 1e7,
                    psi_inf=pinf, ws=pinf / (VT * self.ll), alpha_inf=alpha_inf, Weff_nm=pr['Weff'] * 1e7,
                    profile=dict(x_nm=(x * 1e7).tolist()[::40], psi=ps.tolist()[::40]))


def dL_ward(dev, vg, n=6000):
    """Ward's original Sec. 6.4 evaluation (p. 121 oxide-fringe profile, decay 2 t_ox/pi) for this device."""
    psig, pinf, Wd = dev.far_field(vg)
    dec = 2 * dev.tox / math.pi
    xg = np.concatenate([[0.0], np.geomspace(1e-10, 60 * max(dev.epsL, dec), n)])
    pg = pinf + (dev.psibi - pinf) * np.exp(-xg / dec)
    qinf = dev.qward(pinf, dev.ainf)
    qg = np.array([dev.qward(p, dev.a0(t / dev.epsL)) for p, t in zip(pg, xg)])
    return float(np.trapezoid(1 - qinf / qg, xg)) * 1e7


def ns_1d(dev, psi_s):
    """Electron sheet density (cm^-2) of the semi-infinite 1-D slice at surface potential psi_s."""
    def E(p):
        b = (p - dev.psib) / VT
        return math.sqrt(max(2 * Q * VT * (dev.p0 * (math.exp(-b) + b - 1) + dev.n0 * (math.expm1(b) - b)) / ESI, 1e-300))
    pts = np.linspace(dev.psib + 1e-6, psi_s, 40)
    return sum(quad(lambda p: NI * math.exp(p / VT) / E(p), a, b, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))


def vg_for_ns(dev, ns_target):
    def f(vg):
        _, p, _ = dev.far_field(vg)
        return math.log(ns_1d(dev, p) / ns_target)
    return brentq(f, 0.2, 6.0, xtol=1e-6)


if __name__ == '__main__':
    d = Device()
    for vg in (1.0, 1.25, 1.5, 2.0, 3.0):
        r = d.dL(vg)
        print(vg, {k: round(v, 4) for k, v in r.items() if k != 'profile'})
