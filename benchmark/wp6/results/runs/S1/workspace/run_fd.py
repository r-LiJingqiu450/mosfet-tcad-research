import time, sys
import mosfet_deck as M
import devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    r={}
    for c in ('source','drain','body'):
        r[c]=ds.get_contact_current(device='mos',contact=c,equation='ElectronContinuityEquation')
    r['V1']=ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return r
print(0, cur(), flush=True)
for dv in [1e-4,-1e-4,1e-3,-1e-3]:
    ds.circuit_alter(name='V1', value=dv); h=M.solve_dc()
    print(dv, h['iterations'], cur(), flush=True)
