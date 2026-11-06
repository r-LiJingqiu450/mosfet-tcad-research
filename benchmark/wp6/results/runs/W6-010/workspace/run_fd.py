import time, sys, json
import devsim as ds
import mosfet_deck as M
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    return dict(Ie_drain=ds.get_contact_current(device='mos',contact='drain',equation='ElectronContinuityEquation'),
                Ie_source=ds.get_contact_current(device='mos',contact='source',equation='ElectronContinuityEquation'),
                Ie_body=ds.get_contact_current(device='mos',contact='body',equation='ElectronContinuityEquation'),
                IV1=ds.get_circuit_node_value(solution='dcop', node='V1.I'))
res={}
res[0.0]=cur(); print(0.0,res[0.0],flush=True)
for dv in [1e-4,-1e-4,1e-3,-1e-3]:
    ds.circuit_alter(name='V1', value=dv)
    r=M.solve_dc(); res[dv]=cur(); print(dv,r['iterations'],res[dv],flush=True)
for d in [1e-4,1e-3]:
    for k in res[0.0]:
        print(d,k,(res[d][k]-res[-d][k])/(2*d), (res[d][k]-res[0.0][k])/d)
