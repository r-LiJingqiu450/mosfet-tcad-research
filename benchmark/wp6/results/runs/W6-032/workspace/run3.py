import numpy as np, json
import mosfet_deck as M, devsim as ds
info=M.build()
X=np.array(M.MESH['x_cm'])
xm=0.5*(X[1:]+X[:-1])
ds.edge_from_node_model(device='mos',region='si',node_model='x')
ds.edge_from_node_model(device='mos',region='si',node_model='y')
y0=np.array(ds.get_edge_model_values(device='mos',region='si',name='y@n0'))
y1=np.array(ds.get_edge_model_values(device='mos',region='si',name='y@n1'))
x0=np.array(ds.get_edge_model_values(device='mos',region='si',name='x@n0'))
x1=np.array(ds.get_edge_model_values(device='mos',region='si',name='x@n1'))
C=np.array(ds.get_edge_model_values(device='mos',region='si',name='EdgeCouple'))
print('diag couple max', np.abs(C[(x0!=x1)&(y0!=y1)]).max())
hor=(y0==y1)
col=np.searchsorted(X, np.minimum(x0,x1)+1e-18)-1
def prof():
    J=np.array(ds.get_edge_model_values(device='mos',region='si',name='ElectronCurrent'))
    s=np.sign(x1-x0)
    return np.bincount(col[hor], weights=(J*C*s)[hor], minlength=len(xm))
p0=prof()
r={}
for dv in (1e-4,-1e-4):
    ds.circuit_alter(name='V1',value=dv); M.solve_dc(); r[dv]=prof()
G=(r[1e-4]-r[-1e-4])/2e-4
np.savez('prof.npz',xm=xm,p0=p0,G=G)
for i in range(0,len(xm),8): print(f"{xm[i]*1e4:8.4f} {p0[i]: .3e} {G[i]: .4e}")
