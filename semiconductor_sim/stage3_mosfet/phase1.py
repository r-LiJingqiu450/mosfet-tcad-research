"""Canonical reproducible phase-1 experiment (numerical, not experimental validation)."""
import hashlib,json,platform,subprocess,sys,time,csv
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import brentq
from common import constants as C
from common.jsonio import dumps,write_json
from .config import PRODUCTION,TableConfig,production_metadata
from common.grid import make_surface_grid
from .table import SliceTable
from .gca import gca_solve,id_vg
from .postprocess import extract_vt,continuity_check,edge_continuity_check,extract_ss
from .reference import reference_slice,reference_current,square_law_constant_current,charge_budget
PKG=Path(__file__).resolve().parents[1]
ROOT=PKG.parent
OUT=PKG/'results'/'phase1'
CACHE=ROOT/'.cache'


def cached_table(nps,npn,mode,config=None):
    config=config or TableConfig(nps=nps,npn=npn,interpolation=mode)
    if (nps,npn,mode)!=(config.nps,config.npn,config.interpolation):
        raise ValueError('Contradictory table configuration')
    files=[Path(__file__).with_name(n) for n in ['solver.py','table.py','config.py']]
    files += [PKG/'common'/n for n in ['constants.py','grid.py','thomas.py']]
    key=b''.join(f.read_bytes() for f in files)+dumps(config.to_dict(),sort_keys=True).encode()
    key+=f'{np.__version__}/{scipy.__version__}'.encode()
    token=hashlib.sha256(key).hexdigest()[:16]
    CACHE.mkdir(exist_ok=True)
    path=CACHE/f'table_{nps}_{npn}_{mode}_{token}.npz'
    x=make_surface_grid(config.x_max,d0=config.d0,h_max=config.hmax)
    t=SliceTable(np.linspace(config.phi_s_min,config.phi_s_max,nps),
                 np.linspace(0,config.phi_n_max,npn),x,sheet_interpolation=mode)
    if path.exists():
        with np.load(path,allow_pickle=False) as d:
            if json.loads(str(d['config_json']))!=config.to_dict():
                raise ValueError('Cache configuration mismatch')
            t.vg=d['vg'];t.ns=d['ns'];t.qs=d['qs']
            t.build_stats=json.loads(str(d['build_stats_json']))
        return t,dict(cache_hit=True,cache_key=token)
    start=time.perf_counter();t.build()
    np.savez_compressed(path,vg=t.vg,ns=t.ns,qs=t.qs,
        config_json=dumps(config.to_dict()),build_stats_json=dumps(t.build_stats))
    return t,dict(cache_hit=False,cache_key=token,build_seconds=time.perf_counter()-start)


def production_table():
    return cached_table(PRODUCTION.nps,PRODUCTION.npn,PRODUCTION.interpolation,config=PRODUCTION)


def write_csv(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def run_study(verbose=True,vg_max=2.,vd_max=2.):
    if vg_max!=2 or vd_max!=2:raise ValueError('Frozen phase-1 matrix has max biases 2 V')
    OUT.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    def say(s):
        if verbose:print(s,flush=True)
    say('Building/reusing production 400x220 log table')
    t,perf=production_table()
    say('Building/reusing comparison 240x110 linear table')
    coarse,cp=cached_table(240,110,'linear')
    bp=OUT/'baseline_recomputed.json'
    say('Recomputing unmodified recovered v1 in isolated interpreter')
    r=subprocess.run([sys.executable,str(ROOT/'scripts/baseline_study.py')],capture_output=True,text=True,check=True)
    write_json(bp,json.loads(r.stdout))
    baseline=json.loads(bp.read_text())
    vg=np.linspace(0,2,81);ilin=id_vg(t,vg,.05);isat=id_vg(t,vg,1.5)
    vt,ith=extract_vt(vg,ilin);inv=square_law_constant_current(.05,ith)
    threshold=dict(vt_cc=vt,vt_long=C.VT_LONG,square_law_cc=inv['vg'],branch=inv['branch'],
                   criterion_offset_mV=1000*inv['offset'],remaining_offset_mV=1000*(vt-inv['vg']),
                   total_offset_mV=1000*(vt-C.VT_LONG))
    say('Independent first-integral reference: 42 bias points')
    refs=[]
    for v in [.5,.75,1,1.25,1.5,1.75,2]:
        for d in [.025,.05,.2,.75,1.5,2]:
            ref=reference_current(v,d);cur=gca_solve(v,d,t)['ID'];ci=gca_solve(v,d,coarse)['ID']
            refs.append(dict(vg=v,vd=d,reference_A=ref['ID'],revised_A=cur,
               relative_error=cur/ref['ID']-1,coarse_linear_relative_error=ci/ref['ID']-1,
               reference_quad_estimate_A=ref['quadrature_abs_estimate']))
    continuity=[]
    for n in [81,161,321,641,1281]:
        r=gca_solve(2,1.5,t,y=np.linspace(0,C.L_CH,n))
        continuity.append(dict(ny=n,finite_difference=continuity_check(r)[0],edge_consistency=edge_continuity_check(r)))
    src=t.solve_vg(2,0);rs=reference_slice(2,0)
    budgets=dict(table=charge_budget(2,*src),independent_reference=charge_budget(2,rs['phi_s'],rs['ns'],rs['qs']))
    from .legacy_gca import gca_solve as old_gca
    ablation=[]
    for v,d in [(2,.05),(2,1.5),(1.175,2),(1.2,2)]:
        old=old_gca(v,d,coarse,clm=True);off=old_gca(v,d,coarse,clm=False)
        beta=C.MU_N*C.COX*C.W_CH/C.L_CH;ov=max(v-baseline['vt_cc'],0.)
        square=beta*(ov*d-.5*d*d) if d<ov else .5*beta*ov*ov
        ablation.append(dict(vg=v,vd=d,square_A=square,legacy_clm_A=old['ID'],legacy_no_clm_A=off['ID'],
            full_integral_same_table_A=gca_solve(v,d,coarse)['ID'],revised_production_A=gca_solve(v,d,t)['ID'],
            delta_L_fraction=1-old['y_p']/C.L_CH))
    say('Spatial refinement without table interpolation')
    grid=[]
    from .solver import solve_slice,surface_charge,electron_sheet_density
    for h in [.5,.25,.125,.0625]:
        xx=make_surface_grid(1.2e-6,d0=h*1e-9,h_max=h*32e-9)
        def direct(ps):
            rr=solve_slice(ps,0,xx)
            if not rr['converged']:raise RuntimeError(rr['status'])
            q,_=surface_charge(rr,xx);return ps-q/C.COX-2
        ps=brentq(direct,.8,.9,xtol=1e-13);rr=solve_slice(ps,0,xx);ns=electron_sheet_density(rr,xx)
        grid.append(dict(d0_nm=h,hmax_nm=h*32,nodes=len(xx),phi_s=ps,
                    source_ns_relative_error=ns/rs['ns']-1,correction_inf=rr['correction_inf'],backward_error=rr['backward_error']))
    s=dict(version='phase1-1.2',environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()),
        assumptions='classical Boltzmann, constant mobility, ideal reservoirs, no default CLM; not 2D DD',
        table_config=PRODUCTION.to_dict(),production_metadata=production_metadata(),table_build_stats=t.build_stats,
        threshold=threshold,ss_mVdec=extract_ss(vg,ilin),charge_budget=budgets,current_reference=refs,
        continuity=continuity,ablation=ablation,spatial_refinement=grid,
        maximum_reference_relative_error=max(abs(r['relative_error']) for r in refs),
        maximum_coarse_linear_relative_error=max(abs(r['coarse_linear_relative_error']) for r in refs),
        performance=perf,comparison_table_performance=cp,baseline=baseline,elapsed_seconds=time.perf_counter()-start)
    s['acceptance']=dict(independent_current_42_points_below_0p5pct=s['maximum_reference_relative_error']<.005,
        saturation_fd_641_below_1e_minus_3=continuity[3]['finite_difference']<1e-3,
        edge_internal_consistency_below_1e_minus_9=max(r['edge_consistency'] for r in continuity)<1e-9,
        threshold_correct_branch=inv['branch']=='saturation')
    from .attribution import run_attribution
    s['attribution']=run_attribution(s,t,OUT)
    write_json(OUT/'phase1_results.json',s)
    for name,rows in [('independent_reference',refs),('continuity',continuity),('clm_ablation',ablation),('spatial_refinement',grid)]:write_csv(name+'.csv',rows)
    np.savez_compressed(OUT/'transfer_curves.npz',vg=vg,id_lin=ilin,id_sat=isat)
    from .phase1_plots import make_plots
    make_plots(s,t,OUT)
    say(dumps({k:s[k] for k in ['threshold','maximum_reference_relative_error','acceptance']},indent=2))
    if not all(s['acceptance'].values()):raise RuntimeError('Phase-1 acceptance failed: inspect results JSON')
    return s


if __name__=='__main__':run_study()
