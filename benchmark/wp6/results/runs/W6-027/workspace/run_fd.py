import time, json, sys
import mosfet_deck as M
import devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    I=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    Ide=ds.get_contact_current(device='mos', contact='drain', equation='ElectronContinuityEquation')
    Is=ds.get_contact_current(device='mos', contact='source', equation='ElectronContinuityEquation')
    Ib=ds.get_contact_current(device='mos', contact='body', equation='ElectronContinuityEquation')
    return I, Ide, Is, Ib
res={0.0:cur()}; print(0.0, res[0.0], flush=True)
for dv in [1e-4,-1e-4,1e-3,-1e-3]:
    ds.circuit_alter(name='V1', value=dv)
    r=M.solve_dc(); res[dv]=cur(); print(dv, r['iterations'], res[dv], flush=True)
for h in [1e-4,1e-3]:
    for k in range(4):
        print(h,k,(res[h][k]-res[-h][k])/(2*h))
