import time, json, numpy as np
import mosfet_deck as M, devsim as ds
t=time.time(); info=M.build(); print(info, time.time()-t)
def cur():
    out={}
    for c in ('source','drain','body','gate'):
        try:
            out[c]=ds.get_contact_current(device='mos',contact=c,equation='ElectronContinuityEquation')
        except Exception as e: out[c]=None
    out['V1.I']=ds.get_circuit_node_value(node='V1.I',solution='dcop')
    return out
print(0, cur())
res={}
for dv in (1e-4,-1e-4,1e-3,-1e-3,1e-2,-1e-2):
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc(); c=cur(); res[dv]=c; print(dv, r['iterations'], c, flush=True)
for d in (1e-4,1e-3,1e-2):
    for k in ('drain','V1.I'):
        print(d,k,(res[d][k]-res[-d][k])/(2*d))
ds.circuit_alter(name='V1', value=0.0); M.solve_dc(); print(cur())
for f in (1.0, 1e3):
    ds.solve(type='ac', frequency=f)
    print('ac',f, ds.get_circuit_node_value(node='V1.I',solution='ssac_real'), ds.get_circuit_node_value(node='V1.I',solution='ssac_imag'))
