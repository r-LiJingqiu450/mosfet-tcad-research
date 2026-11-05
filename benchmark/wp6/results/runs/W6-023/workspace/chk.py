import numpy as np, devsim as ds, mosfet_deck as M
M.build()
ds.circuit_alter(name='V1', value=1e-3); r=M.solve_dc()
g=lambda r,n: np.array(ds.get_node_model_values(device=M.DEV, region=r, name=n))
for reg in (M.SI, M.OX):
    x,y=M.coords(reg); psi=g(reg,'Potential'); lat=g(reg,'PsiLat')
    for c in ('drain','oxc_d','source','oxc_s'):
        try: idx=ds.get_contact_nodes if False else None
        except: pass
    m=(np.abs(x-x.max())<1e-15)&(y<=M.MESH['contact_extent_cm']*1.000001)&(y>=-1e-15)
    d=psi[m]-lat[m]
    print(reg,'drain-side nodes',m.sum(),'psi-lat min/max',d.min(),d.max(), 'y of extremes', y[m][np.argmin(d)], y[m][np.argmax(d)])
    m=(np.abs(x-x.min())<1e-15)&(y<=M.MESH['contact_extent_cm']*1.000001)&(y>=-1e-15)
    d=psi[m]-lat[m]
    print(reg,'source-side nodes',m.sum(),'psi-lat min/max',d.min(),d.max())
    if reg==M.SI:
        n=g(reg,'Electrons'); m=(np.abs(x-x.max())<1e-15)&(y<=M.MESH['contact_extent_cm']*1.000001)
        print('n/(ni exp((psi-V)/Vt))-1', np.abs(n[m]/(M.NI*np.exp((psi[m]-1e-3)/M.VT))-1).max())
        x0=np.abs(x-x.max())<1e-15; print('Si y=0 corner psi', psi[x0&(np.abs(y)<1e-15)], lat[x0&(np.abs(y)<1e-15)])
    else:
        x0=np.abs(x-x.max())<1e-15; print('Ox corner', psi[x0&(np.abs(y)<1e-15)], lat[x0&(np.abs(y)<1e-15)])
print(ds.get_parameter(device=M.DEV,name='s_couple'))
