import time, sys
import devsim as ds
import mosfet_deck as M
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def I(c):
    return ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation')
def report(V):
    cur={c:I(c) for c in ('source','drain','body')}
    print(V, cur, 'circ', ds.get_circuit_node_value(solution='dcop', node='V1.I'), flush=True)
    return cur
I0=report(0.0)
# AC small-signal
for f in (1e-3, 1.0):
    ds.solve(type='ac', frequency=f)
    print('ac', f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'), flush=True)
res={}
for dv in (1e-4,-1e-4,1e-3,-1e-3):
    ds.circuit_alter(name='V1', value=dv); M.solve_dc(); res[dv]=report(dv)
for h in (1e-4,1e-3):
    print('FD h',h,(res[h]['drain']-res[-h]['drain'])/(2*h))
