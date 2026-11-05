import time, json, numpy as np
import mosfet_deck as M, devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t)
def cur():
    return dict(drain=ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation'),
                source=ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation'),
                body=ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation'),
                circ=ds.get_circuit_node_value(solution='dcop', node='V1.I'))
I0=cur(); print('V=0', I0)
res={}
for dv in [1e-3, 3e-4, 1e-4, 3e-5]:
    r={}
    for s in (+1,-1):
        ds.circuit_alter(name='V1', value=s*dv); h=M.solve_dc(); r[s]=cur(); print(s*dv, h['iterations'], r[s], flush=True)
    for k in r[1]:
        print(' dv',dv,k,'G=',(r[1][k]-r[-1][k])/(2*dv))
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in [1e-3, 1.0, 1e3]:
    ds.solve(type='ac', frequency=f)
    print('AC f',f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
