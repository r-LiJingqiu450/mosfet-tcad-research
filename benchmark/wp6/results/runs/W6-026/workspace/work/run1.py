import sys, time; sys.path.insert(0,'.')
import devsim as ds, mosfet_deck as M
t=time.time(); info=M.build(); print(info, time.time()-t)
def cur():
    out={}
    for c in ('source','drain','body'):
        out[c]=ds.get_contact_current(device='mos',contact=c,equation='ElectronContinuityEquation')
    out['V1']=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return out
print(0, cur())
for dv in (1e-4,-1e-4,1e-3,-1e-3):
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc(); print(dv, r['iterations'], cur())
ds.circuit_alter(name='V1', value=0.0); M.solve_dc(); print(0, cur())
ds.solve(type='ac', frequency=1.0)
print('ac', ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
