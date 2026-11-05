import sys, time, json
import numpy as np
import devsim as ds
import mosfet_deck as M
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    Ic=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    Id=ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation')
    Is=ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation')
    Ib=ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation')
    return Ic,Id,Is,Ib
print('V=0', cur(), flush=True)
res={}
for dv in [1e-3,3e-4,1e-4,-1e-4,-3e-4,-1e-3]:
    ds.circuit_alter(name='V1', value=dv)
    r=M.solve_dc(); c=cur(); res[dv]=c
    print(dv, r['iterations'], c, flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
I0=cur(); print('V=0 again', I0)
for dv in [1e-3,3e-4,1e-4]:
    print('central dv',dv, [(res[dv][k]-res[-dv][k])/(2*dv) for k in range(4)])
    print('fwd dv',dv, [(res[dv][k]-I0[k])/(dv) for k in range(4)])
