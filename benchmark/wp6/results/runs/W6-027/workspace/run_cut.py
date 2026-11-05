import time, json, sys
import numpy as np
import mosfet_deck as M
import devsim as ds
D,R='mos','si'
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
ds.edge_from_node_model(device=D, region=R, node_model='node_index')
i0=np.array(ds.get_edge_model_values(device=D,region=R,name='node_index@n0')).astype(int)
i1=np.array(ds.get_edge_model_values(device=D,region=R,name='node_index@n1')).astype(int)
C=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
x,y=M.coords(R)
Lmid=1e-4; ycut=1.1e-4
def wfun(xc, yc):
    return ((x>xc)&(y<yc)).astype(float)
cuts=[(1.0e-4,1.1e-4),(0.8e-4,1.0e-4),(1.3e-4,1.15e-4)]
W=[wfun(*c) for c in cuts]
def cur():
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))
    F=[float(np.sum((w[i0]-w[i1])*J*C)) for w in W]
    I=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    Ide=ds.get_contact_current(device=D, contact='drain', equation='ElectronContinuityEquation')
    return [I,Ide]+F
res={}
for dv in [0.0,1e-4,-1e-4,1e-3,-1e-3,2e-3,-2e-3,0.05]:
    ds.circuit_alter(name='V1', value=dv)
    r=M.solve_dc(); res[dv]=cur(); print(dv, r['iterations'], res[dv], flush=True)
for h in [1e-4,1e-3,2e-3]:
    print('h',h,[(a-b)/(2*h) for a,b in zip(res[h],res[-h])])
np.save('res_cut.npy', np.array([[k]+v for k,v in res.items()]))
