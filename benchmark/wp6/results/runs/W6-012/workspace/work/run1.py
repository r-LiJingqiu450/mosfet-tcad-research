import sys, time; sys.path.insert(0,'.')
import mosfet_deck as M, devsim as ds
t=time.time(); info=M.build(); print(info, time.time()-t)
def cur():
    return {c: ds.get_contact_current(device='mos',contact=c,equation='ElectronContinuityEquation') for c in ('source','drain','body')}, ds.get_circuit_node_value(solution='dcop', node='V1.I')
print(0, cur())
for v in (1e-4,-1e-4,1e-3,-1e-3):
    ds.circuit_alter(name='V1', value=v); r=M.solve_dc(); print(v, r['iterations'], cur())
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in (1e-3,1.0):
    ds.solve(type='ac', frequency=f)
    print('ac',f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
