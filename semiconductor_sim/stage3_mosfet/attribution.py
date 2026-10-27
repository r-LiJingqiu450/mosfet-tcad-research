"""Paired numerical attribution. Signed contributions add to total error.

Direct finite-difference gate solves remove table interpolation. Reference-based
surrogate tables separate phi_s/VG inversion interpolation from phi_n interpolation.
These are conditional numerical ablations, not unique physical decompositions.
"""
from functools import lru_cache
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from common import constants as C
from common.grid import make_surface_grid
from common.jsonio import write_json
from .solver import solve_poisson_slice,surface_charge,electron_sheet_density
from .reference import reference_slice,reference_current
from .gca import gca_solve,_areas


class DirectSlices:
    def __init__(self,d0=.125e-9,hmax=4e-9,clip=None):
        self.x=make_surface_grid(1.2e-6,d0=d0,h_max=hmax)
        self.clip=clip
        self.minimum=float('inf');self.maximum=float('-inf');self.clip_evaluations=0

    @lru_cache(maxsize=8192)
    def sheet(self,vg,pn):
        def state(ps):
            r=solve_poisson_slice(ps,pn,self.x,density_clip_limit=self.clip)
            if not r['converged']:raise RuntimeError(r['status'])
            self.minimum=min(self.minimum,r['potential_min'])
            self.maximum=max(self.maximum,r['potential_max'])
            self.clip_evaluations+=r['clip_evaluations']
            return r
        def gate(ps):
            q,_=surface_charge(state(ps),self.x)
            return C.VFB+ps-q/C.COX-vg
        ps=brentq(gate,0.,min(.9+pn,2.05),xtol=2e-13)
        return electron_sheet_density(state(ps),self.x)

    def current(self,vg,vd,epsrel=2e-7):
        val,err=quad(lambda pn:self.sheet(float(vg),float(pn)),0,vd,
                     epsabs=1e-4,epsrel=epsrel,limit=200)
        pref=C.Q_E*C.W_CH*C.MU_N/C.L_CH
        return float(pref*val),float(pref*err)


def run_attribution(study,table,out):
    worst=max(study['current_reference'],key=lambda r:abs(r['relative_error']))
    vg,vd=worst['vg'],worst['vd'];ref=reference_current(vg,vd)
    direct=DirectSlices();idir,quaderr=direct.current(vg,vd)
    # Same grid and algorithm, historical clamp on/off, cached separately.
    clipped=DirectSlices(clip=3.);iclip,_=clipped.current(vg,vd)
    iref=ref['ID'];iprod=gca_solve(vg,vd,table)['ID']
    # Exact interpolation of the independent continuum Ns at production phi_n knots.
    # This estimates phi_n-only interpolation on the continuum reference.
    pn=np.r_[0.,table.phi_n_grid[(table.phi_n_grid>0)&(table.phi_n_grid<vd)],vd]
    ns=np.array([reference_slice(vg,float(v))['ns'] for v in pn])
    ipn=C.Q_E*C.W_CH*C.MU_N/C.L_CH*float(np.sum(_areas(pn,ns,'log')))
    grid=(idir-iref)/iref
    interp=(iprod-idir)/iref
    clip=(iclip-idir)/iref
    # Production is unclipped. The clamp contribution is a paired sensitivity,
    # not added into the production decomposition a second time.
    attribution=dict(point=dict(vg=vg,vd=vd),reference_A=iref,
        production_A=iprod,direct_same_grid_A=idir,clipped_direct_A=iclip,
        total_relative=(iprod-iref)/iref,grid_and_boundary_relative=grid,
        table_interpolation_relative=interp,
        telescoping_closure_relative=(iprod-iref)/iref-grid-interp,
        clamp_paired_relative=clip,clamp_activated_values=clipped.clip_evaluations,
        potential_min_V=direct.minimum,potential_max_V=direct.maximum,
        phi_n_only_interpolation_relative=(ipn-iref)/iref,
        direct_quad_estimate_relative=quaderr/iref,
        evidence_type='conditional_numerical_ablation')
    # Refinement checks whether the grid+boundary component behaves as discretization.
    levels=[]
    for nm in [.25,.125,.0625]:
        obj=direct if nm==.125 else DirectSlices(d0=nm*1e-9,hmax=nm*32e-9)
        cur,err=obj.current(vg,vd)
        levels.append(dict(d0_nm=nm,hmax_nm=32*nm,current_A=cur,
                           relative_to_reference=cur/iref-1,quad_estimate_A=err))
    backgrounds=[]
    for gate in [.05,.1,.25,.5,.75,1.]:
        for drain in [.025,.05,1.5]:
            r=reference_current(gate,drain)
            backgrounds.append(dict(vg=gate,vd=drain,total_A=r['ID'],
                 excess_A=r['ID_excess'],neutral_background_A=r['ID_background'],
                 neutral_background_fraction=r['background_fraction']))
    report=dict(worst_point=attribution,direct_grid_refinement=levels,
                neutral_background=backgrounds,
                table_trial_potential_range=study['table_build_stats'],
                interpretation='background is a bookkeeping component in the reduced model; not experimental current validation')
    write_json(out/'error_attribution.json',report)
    import csv
    with (out/'neutral_background.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(backgrounds[0]));w.writeheader();w.writerows(backgrounds)
    with (out/'error_attribution.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=['component','relative_contribution'])
        w.writeheader()
        for key in ['total_relative','grid_and_boundary_relative','table_interpolation_relative',
                    'clamp_paired_relative','phi_n_only_interpolation_relative','direct_quad_estimate_relative']:
            w.writerow(dict(component=key,relative_contribution=attribution[key]))
    return report
