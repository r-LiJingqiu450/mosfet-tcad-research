"""Deterministic scientific plots; SVG and PNG outputs, English labels."""
import numpy as np
from .gca import id_vg


def make_plots(s,t,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25})
    def save(f,name):
        f.tight_layout();f.savefig(out/(name+'.png'),dpi=180);f.savefig(out/(name+'.svg'));plt.close(f)
    old=s['baseline']['continuity'];new=s['continuity'];f,ax=plt.subplots(figsize=(7,4.5))
    ax.loglog([r['ny'] for r in old],[r['finite_difference'] for r in old],'o-',label='Recovered v1 + CLM')
    ax.loglog([r['ny'] for r in new],[r['finite_difference'] for r in new],'s-',label='Revised full integral (FD check)')
    ax.axhline(1e-3,color='black',ls='--',label='Original acceptance target')
    ax.set(xlabel='Output y nodes',ylabel='Relative point-current fluctuation');ax.legend();save(f,'fig1_continuity')
    rows=s['current_reference'];vs=sorted({r['vg'] for r in rows});ds=sorted({r['vd'] for r in rows})
    f,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,key,title in zip(axes,['coarse_linear_relative_error','relative_error'],['240x110 linear table','400x220 log table']):
        vals=np.array([[next(r[key] for r in rows if r['vg']==v and r['vd']==d)*100 for d in ds] for v in vs])
        im=ax.imshow(vals,origin='lower',aspect='auto',vmin=-.5,vmax=5,cmap='coolwarm')
        ax.set(xticks=range(len(ds)),xticklabels=ds,yticks=range(len(vs)),yticklabels=vs,xlabel='VD (V)',ylabel='VG (V)',title=title)
        ax.grid(False);f.colorbar(im,ax=ax,label='Current error (%)')
    save(f,'fig2_independent_reference')
    b=s['ablation'][1];f,ax=plt.subplots(figsize=(7,4.5))
    vals=[b['square_A'],b['legacy_no_clm_A'],b['legacy_clm_A'],b['revised_production_A']]
    ax.bar(['Square law','Legacy no CLM','Legacy + CLM','Revised no CLM'],np.array(vals)*1e6)
    ax.set(ylabel='ID (uA)',title='VG=2 V, VD=1.5 V: CLM ablation');save(f,'fig3_clm_ablation')
    data=np.load(out/'transfer_curves.npz');f,axes=plt.subplots(1,2,figsize=(10,4));base=s['baseline']
    axes[0].semilogy(base['vg'],base['transfer_legacy'],label='Recovered v1 + CLM')
    axes[0].semilogy(data['vg'],data['id_sat'],label='Revised classical, no CLM')
    axes[0].set(xlabel='VG (V)',ylabel='ID (A)',ylim=(1e-16,2e-4));axes[0].legend()
    v=np.linspace(1.1,1.3,201);axes[1].plot(v,id_vg(t,v,2)*1e6)
    axes[1].set(xlabel='VG (V)',ylabel='ID (uA)',title='No hard CLM switch near VG=1.2 V');save(f,'fig4_transfer')
    gr=s['spatial_refinement'];f,ax=plt.subplots(figsize=(7,4.5))
    ax.loglog([r['d0_nm'] for r in gr],[abs(r['source_ns_relative_error']) for r in gr],'o-')
    ax.set(xlabel='Surface spacing (nm); hmax=32 h0',ylabel='Source Ns error vs independent reference',title='Spatial refinement only, VG=2 V');save(f,'fig5_spatial_refinement')

    # These two figures expose the physical bookkeeping and signed numerical ablations.
    bg=[r for r in s['attribution']['neutral_background'] if r['vd']==.05]
    f,ax=plt.subplots(figsize=(7,4.5))
    ax.semilogy([r['vg'] for r in bg],[100*r['neutral_background_fraction'] for r in bg],'o-')
    ax.set(xlabel='VG (V)',ylabel='Neutral-background fraction (%)',
           title='Reduced-model background contribution, VD=0.05 V')
    save(f,'fig6_background')
    a=s['attribution']['worst_point'];f,ax=plt.subplots(figsize=(8,4.5))
    labels=['Grid + boundary','Table interpolation','Total error','Clamp sensitivity']
    vals=[a[k]*100 for k in ['grid_and_boundary_relative','table_interpolation_relative','total_relative','clamp_paired_relative']]
    bars=ax.bar(labels,vals);ax.axhline(0,color='black',lw=.8)
    ax.bar_label(bars,fmt='%.4f',padding=3)
    ax.set(ylabel='Signed relative contribution (%)',
           title=f"Worst point: VG={a['point']['vg']} V, VD={a['point']['vd']} V")
    save(f,'fig7_error_attribution')
