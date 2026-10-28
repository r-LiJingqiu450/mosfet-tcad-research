#!/usr/bin/env python3
"""Independent S2 Gummel/Poisson convergence trace.

This script mirrors the diagnostic path without changing production solver
code.  It writes the per-outer and per-inner iteration evidence consumed by
the Phase 2A audit report.
"""

import json, csv, warnings, sys, hashlib, platform
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'semiconductor_sim'))
from common.constants import EPS_SI, NI, Q_E, VT, NA_CM3, ND_CM3
from common.grid import make_junction_grid
from stage2_pn import solver as S

Na=NA_CM3*1e6; Nd=ND_CM3*1e6
x=make_junction_grid(2.0e-6)
Na_arr=np.full_like(x,Na); Nd_arr=np.full_like(x,Nd)
Na_arr[x>=0]=0.0; Nd_arr[x<0]=0.0
bc=S.contact_values(0.0,Na,Nd)
cc=dict(Na=Na,Nd=Nd)

# Exact mirror of solve_poisson_pn, with per-iteration diagnostics.
def poisson_diag(phi,n,p,tol=1e-11,max_iter=100,clamp=0.3,damped=True,max_backtrack=12):
    hw=np.diff(x); hL=hw[:-1]; hR=hw[1:]; h_i=0.5*(hL+hR)
    aa=-EPS_SI/hL; bb0=EPS_SI/hL+EPS_SI/hR; cc0=-EPS_SI/hR
    ph=phi.astype(float).copy(); ph[0]=bc['phi_left']; ph[-1]=bc['phi_right']
    F=S._poisson_residual(ph,n,p,Na_arr,Nd_arr,aa,bb0,cc0,h_i)
    normF=float(np.linalg.norm(F)); initial=normF; rows=[]; armijo_fails=0
    converged=False; final_alpha=1.0; final_dmove=float('nan')
    for it in range(1,max_iter+1):
        b=bb0.copy(); a=aa.copy(); c=cc0.copy(); a[0]=0.; c[-1]=0.
        delta=S.thomas_solve(a,b,c,-F); delta=np.clip(delta,-clamp,clamp)
        delta_inf=float(np.max(np.abs(delta))); alpha=1.0; norm_try=normF; backtracks=0; accepted=False
        for j in range(max(1,int(max_backtrack))):
            ph_try=ph.copy(); ph_try[1:-1]+=alpha*delta
            F_try=S._poisson_residual(ph_try,n,p,Na_arr,Nd_arr,aa,bb0,cc0,h_i)
            norm_try=float(np.linalg.norm(F_try)); backtracks=j
            if norm_try <= (1.0-2e-4*alpha)*normF:
                accepted=True; break
            alpha*=0.5
        if not accepted: armijo_fails += 1
        src_scale=Q_E*max(np.max(Na_arr),np.max(Nd_arr))*np.mean(h_i)*np.sqrt(max(len(ph)-2,1))
        floor_exempt = bool(normF <= 1e-8*src_scale)
        ph,F,norm_prev=ph_try,F_try,normF
        normF=norm_try
        dmove=float(alpha*delta_inf)
        rows.append(dict(kind='newton',iteration=it,norm_before=float(norm_prev),norm_after=normF,relative=normF/max(initial,1e-300),delta_inf=delta_inf,alpha=float(alpha),backtracks=int(backtracks+1),armijo_accepted=accepted,source_scale=float(src_scale),floor_exempt=floor_exempt,dmove=dmove))
        final_alpha=alpha; final_dmove=dmove
        if dmove<tol:
            converged=True; break
    return ph,int(it),bool(converged),rows,dict(final_norm=float(normF),initial_norm=float(initial),armijo_failures=int(armijo_fails),final_alpha=float(final_alpha),final_dmove=float(final_dmove))

def gummel_diag(omega,max_outer=30,tol=1e-9):
    phi=S.initial_equilibrium_guess(x,Na,Nd)
    phi=phi+(bc['phi_left']-(-VT*np.log(Na/NI)))
    phi=phi.astype(float).copy(); phi[0]=bc['phi_left']; phi[-1]=bc['phi_right']
    ss=(x-x[0])/(x[-1]-x[0]); phin=bc['phi_n_L']*(1-ss)+bc['phi_n_R']*ss; phip=bc['phi_p_L']*(1-ss)+bc['phi_p_R']*ss
    n,p=S.carrier_from_phi(phi,phin,phip); n=n.astype(float).copy(); p=p.astype(float).copy(); n[0],p[0],n[-1],p[-1]=S.contact_carriers(bc,Na,Nd)
    outer=[]; all_newton=[]; pk=dict(tol=1e-11,max_iter=100,clamp=0.3,damped=True)
    for out in range(1,max_outer+1):
        n_s=S.solve_electron_continuity(phi,p,n,x,bc,cc); n_new=(1-omega)*n+omega*n_s
        p_s=S.solve_hole_continuity(phi,n_new,p,x,bc,cc); p_new=(1-omega)*p+omega*p_s
        phi_s,ni,inner_ok, hist, summ=poisson_diag(phi,n_new,p_new,**pk)
        phi_new=(1-omega)*phi+omega*phi_s
        dphi=float(np.max(np.abs(phi_new-phi))); dn=float(np.max(np.abs(n_new-n))/max(np.max(n_new),1.0)); dp=float(np.max(np.abs(p_new-p))/max(np.max(p_new),1.0))
        outer.append(dict(kind='gummel',omega=omega,iteration=out,dphi=dphi,dn=dn,dp=dp,inner_iters=ni,inner_converged=inner_ok,inner_final_norm=summ['final_norm'],inner_initial_norm=summ['initial_norm'],inner_armijo_failures=summ['armijo_failures'],inner_final_alpha=summ['final_alpha'],inner_final_dmove=summ['final_dmove']))
        for h in hist:
            h.update(omega=omega,outer_iteration=out)
        all_newton.extend(hist)
        phi,n,p=phi_new,n_new,p_new
    return outer,all_newton

rows=[]
for om in (1.0,0.5):
    o,n=gummel_diag(om)
    rows += o+n
payload={'metadata':{
    'command':'MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/diagnose_s2_convergence.py',
    'python':sys.version,
    'platform':platform.platform(),
    'solver_sha256':hashlib.sha256((ROOT / 'semiconductor_sim/stage2_pn/solver.py').read_bytes()).hexdigest(),
    'production_solver_modified':False,
},'grid':{'N':len(x),'L_m':2e-6,'min_dx_m':float(np.min(np.diff(x))),'max_dx_m':float(np.max(np.diff(x)))},'params':{'Na_m3':Na,'Nd_m3':Nd,'omega':[1.0,0.5],'outer_max':30,'outer_tol':1e-9,'inner_tol':1e-11,'inner_max_iter':100,'inner_clamp_V':0.3,'armijo_backtrack_max':12},'runs':{}}
for om in (1.0,0.5):
 o=[r for r in rows if r.get('kind')=='gummel' and r['omega']==om]
 nn=[r for r in rows if r.get('kind')=='newton' and r['omega']==om]
 sample=[r for r in nn if r['outer_iteration'] in (1,2,30)
         and r['iteration'] in (1,2,3,98,99,100)]
 payload['runs'][str(om)]={'outer':o,'newton_sample':sample,
                           'summary':{'outer_final_dphi':o[-1]['dphi'],
                                      'outer_final_dn':o[-1]['dn'],
                                      'outer_final_dp':o[-1]['dp'],
                                      'outer_converged':any(r['dphi']<1e-9 and r['iteration']>1 for r in o),
                                      'newton_trace_steps':len(nn),
                                      'newton_nonconverged_calls':sum(not r['inner_converged'] for r in o),
                                      'armijo_failed_steps':sum(r['inner_armijo_failures'] for r in o)}}
(ROOT / 'verification' / 's2_convergence_diagnostics.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
with (ROOT / 'verification' / 's2_convergence_history.csv').open('w', newline='') as f:
    fields=['kind','omega','outer_iteration','iteration','dphi','dn','dp','inner_iters','inner_converged','inner_initial_norm','inner_final_norm','inner_armijo_failures','inner_final_alpha','inner_final_dmove','norm_before','norm_after','relative','delta_inf','alpha','backtracks','armijo_accepted','source_scale','floor_exempt','dmove']
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
    for r in rows:w.writerow({k:r.get(k,'') for k in fields})
for om in (1.0,0.5):
    o=payload['runs'][str(om)]['outer']; print('OMEGA',om,'outer_final',o[-1]); print('first_outer',o[:3]); print('last_outer',o[-5:]);
    nn=[r for r in rows if r.get('kind')=='newton' and r['omega']==om]
    bad=[r for r in nn if not r['armijo_accepted']]
    print('newton steps',len(nn),'armijo_failed',len(bad),'nonconv calls',payload['runs'][str(om)]['summary']['newton_nonconverged_calls']); print('bad sample',bad[:3]);
