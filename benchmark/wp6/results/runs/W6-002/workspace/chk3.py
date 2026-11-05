import numpy as np, devsim as ds, mosfet_deck as M
M.build()
def g(r,n): return np.array(ds.get_node_model_values(device=M.DEV, region=r, name=n))
xs,ys=M.coords(M.SI); xo,yo=M.coords(M.OX)
P0=g(M.SI,'Potential'); Po0=g(M.OX,'Potential'); N0=g(M.SI,'Electrons'); L0=g(M.SI,'PsiLat'); Lo=g(M.OX,'PsiLat')
d=1e-3
ds.circuit_alter(name='V1', value=d); M.solve_dc()
P=g(M.SI,'Potential'); Po=g(M.OX,'Potential'); N=g(M.SI,'Electrons')
dr=(np.abs(xs-xs.max())<1e-15)&(ys<=M.MESH['contact_extent_cm']*(1+1e-12))
sr=(np.abs(xs-xs.min())<1e-15)&(ys<=M.MESH['contact_extent_cm']*(1+1e-12))
print('drain nodes', dr.sum(), 'max|psi-(lat+d)|', np.abs(P[dr]-L0[dr]-d).max(), 'n rel', np.abs(N[dr]/(M.NI*np.exp((P[dr]-d)/M.VT))-1).max())
print('source', np.abs(P[sr]-L0[sr]).max())
for cx in (xs.min(), xs.max()):
    i=np.where((np.abs(xs-cx)<1e-15)&(np.abs(ys)<1e-15))[0]; j=np.where((np.abs(xo-cx)<1e-15)&(np.abs(yo)<1e-15))[0]
    print('corner x',cx,'si psi',P[i],'lat',L0[i],'ox psi',Po[j],'ox lat',Lo[j], 'si0',P0[i],'ox0',Po0[j])
# drain-side well potential change
print('dpsi near drain top', (P-P0)[dr][:5])
