"""Independent classical uniform-substrate Poisson first-integral reference.

Uses neither the finite-difference Poisson solver nor SliceTable. A semi-infinite
neutral bulk is used; finite-domain Ns includes neutral bulk background over
x_max. Its body potential differs from imposed zero by <= ni^2/NA^2 * UT for
nonzero phi_n (about 3e-14 V at default doping). This approximation and finite
bulk treatment are recorded in the report; this is code verification, not
experimental validation. Only positive surface band bending is supported.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from common import constants as C


def _em1mx(z):
    if abs(z) < 1e-3:
        return z*z*(.5+z*(1/6+z*(1/24+z*(1/120+z/720))))
    return np.expm1(z)-z


def reference_slice(vg, phi_n, na=C.NA, x_max=1.2e-6):
    if vg < .05 or phi_n < 0:
        raise ValueError('Reference supports positive inversion/depletion gate bias')
    p0=.5*(na+np.sqrt(na*na+4*C.NI**2)); n0=C.NI**2/p0
    product=C.NI**2*np.exp(-phi_n/C.VT)
    pb=.5*(na+np.sqrt(na*na+4*product)); nb=product/pb
    bulk=-C.VT*np.log(pb/p0)
    def energy(u):
        return C.VT*(pb*_em1mx(-u)+nb*_em1mx(u))
    def charge(ps):
        u=(ps-bulk)/C.VT
        return -np.sqrt(max(0.,2*C.EPS_SI*C.Q_E*energy(u)))
    ps=brentq(lambda v:C.VFB+v-charge(v)/C.COX-vg,
              bulk+1e-14,max(2.5,phi_n+1),xtol=2e-14)
    us=(ps-bulk)/C.VT
    def sheet_integrand(u):
        if abs(u)<1e-8:
            return nb*np.sqrt(C.EPS_SI*C.VT/(C.Q_E*(pb+nb)))
        return nb*np.expm1(u)*C.VT/np.sqrt(2*C.Q_E/C.EPS_SI*energy(u))
    excess,error=quad(sheet_integrand,0,us,epsabs=1e-6,epsrel=2e-10,limit=200)
    return dict(phi_s=ps,qs=charge(ps),ns=excess+nb*x_max,ns_excess=excess,
                ns_background=nb*x_max,nb=nb,
                quadrature_abs_estimate=error,bulk_shift=bulk)


def reference_current(vg, vd, x_max=1.2e-6, epsrel=2e-8):
    excess,error=quad(lambda pn:reference_slice(vg,pn,x_max=x_max)['ns_excess'],0,vd,
                      epsabs=1e-4,epsrel=epsrel,limit=160)
    # Integrate the neutral background independently so its contribution is visible.
    background,bg_error=quad(lambda pn:reference_slice(vg,pn,x_max=x_max)['ns_background'],
                            0,vd,epsabs=1e-6,epsrel=epsrel,limit=160)
    pref=C.Q_E*C.W_CH*C.MU_N/C.L_CH
    return dict(ID=pref*(excess+background),ID_excess=pref*excess,
                ID_background=pref*background,background_fraction=background/(excess+background),
                quadrature_abs_estimate=pref*(error+bg_error))


def square_law_constant_current(vd, ith=1e-7, vt=C.VT_LONG):
    """Exact inverse of the same piecewise square law used by constants.py."""
    if vd <= 0 or ith <= 0:
        raise ValueError('VD and threshold current must be positive')
    beta=C.MU_N*C.COX*C.W_CH/C.L_CH
    if ith <= .5*beta*vd*vd:
        ov=np.sqrt(2*ith/beta);branch='saturation'
    else:
        ov=vd/2+ith/(beta*vd);branch='linear'
    return dict(vg=vt+ov,offset=ov,branch=branch,beta=beta)


def charge_budget(vg, phi_s, ns, qs):
    """Algebraic comparison to a specified square-law reference, not validation.

    Since Qb is defined from Qs and Ns, closure follows from the gate equation;
    its small roundoff error is not evidence of physical accuracy.
    """
    qb_abs=-(qs+C.Q_E*ns)
    qbth=C.COX*(C.VT_LONG-C.VFB-C.TWICE_PSI_F)
    surface=phi_s-C.TWICE_PSI_F
    body=(qb_abs-qbth)/C.COX
    total=vg-C.VT_LONG-C.Q_E*ns/C.COX
    return dict(evidence_type="algebraic_reference_comparison",surface_term_V=surface,body_term_V=body,total_deficit_V=total,
                closure_V=total-surface-body,qb_abs=qb_abs)
