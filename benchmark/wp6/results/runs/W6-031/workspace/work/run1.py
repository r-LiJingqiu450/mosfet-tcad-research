import sys, time, json
sys.path.insert(0, '/home/claude/wp6_runs/W6-031')
import numpy as np, devsim as ds
import mosfet_deck as M
t=time.time(); info=M.build(); print(info, time.time()-t, flush=True)
def currents():
    out={}
    for c in ('source','drain','body'):
        out[c]=ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation')
    out['V1.I']=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return out
print('V=0', currents(), flush=True)
# AC small-signal
for f in (1e-3, 1.0, 1e3):
    ds.solve(type='ac', frequency=f)
    print('AC f', f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'), flush=True)
res={}
for v in (1e-4, -1e-4, 1e-3, -1e-3, 4e-3, -4e-3):
    ds.circuit_alter(name='V1', value=v); M.solve_dc()
    res[v]=currents(); print(v, res[v], flush=True)
for h in (1e-4,1e-3,4e-3):
    for k in ('drain','V1.I','source','body'):
        print(h, k, (res[h][k]-res[-h][k])/(2*h))
