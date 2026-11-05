import time, json
import mosfet_deck as M
import devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    r={}
    for c in ('drain','source','body'):
        r[c]=ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation')
    r['V1.I']=ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return r
I0=cur(); print('V=0', I0, flush=True)
# AC small-signal at low frequency
for f in (1.0,):
    ds.solve(type='ac', frequency=f)
    print('ac', f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'), flush=True)
res={}
for dv in (1e-3, 1e-4, 1e-5):
    out=[]
    for s in (1,-1):
        ds.circuit_alter(name='V1', value=s*dv)
        h=M.solve_dc()
        out.append(cur())
        print(s*dv, len(h['history']), out[-1], flush=True)
    g={k:(out[0][k]-out[1][k])/(2*dv) for k in out[0]}
    print('G', dv, g, flush=True)
    res[dv]=g
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
print('back', cur())
