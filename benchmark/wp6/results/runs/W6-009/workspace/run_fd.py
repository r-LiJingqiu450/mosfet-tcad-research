import time, json
import numpy as np
import mosfet_deck as M
import devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t)
def cur():
    Ic=ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation')
    Is=ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation')
    Ib=ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation')
    Iv=ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return Ic,Is,Ib,Iv
print('V=0', cur())
res={}
for dv in [1e-3,-1e-3,1e-4,-1e-4,1e-5,-1e-5]:
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc(); res[dv]=cur(); print(dv, r['iterations'], res[dv])
for d in [1e-3,1e-4,1e-5]:
    print(d, [(res[d][k]-res[-d][k])/(2*d) for k in range(4)])
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
print('V=0 again', cur())
for f in [0.0, 1.0]:
    try:
        ds.solve(type='ac', frequency=f)
        print('ac', f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'))
    except Exception as e: print('ac fail', f, e)
