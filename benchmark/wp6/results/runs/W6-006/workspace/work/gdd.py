import sys, time, json
sys.path.insert(0, '.')
import mosfet_deck as M
import devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def currents():
    out={}
    for c in ('drain','source','body'):
        out[c]=ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation')
    out['V1.I']=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return out
I0=currents(); print('V=0', I0, flush=True)
res={'0':I0}
for dv in (1e-4,-1e-4,1e-3,-1e-3):
    ds.circuit_alter(name='V1', value=dv)
    r=M.solve_dc(); I=currents(); res[str(dv)]=I
    print(dv, r['iterations'], I, flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for dv in ('0.0001','0.001'):
    for k in I0:
        print(dv,k,'G_central=', (res[dv][k]-res['-'+dv][k])/(2*float(dv)))
# AC small-signal
try:
    ds.solve(type='ac', frequency=1.0)
    print('AC V1.I', ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
    for c in ('drain','source','body'):
        print('AC', c, ds.get_contact_charge  and None)
except Exception as e: print('AC failed', e)
json.dump(res, open('work/fd.json','w'), indent=1)
