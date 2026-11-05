import sys, time, json
import mosfet_deck as M, devsim as ds
t=time.time(); info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    d={}
    for c in ('source','drain','body'):
        d[c]=ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation')
    d['V1.I']=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return d
I0=cur(); print('V=0', I0, flush=True)
res={}
for dv in [float(a) for a in sys.argv[1:]]:
    out={}
    for s in (+1,-1):
        ds.circuit_alter(name='V1', value=s*dv); r=M.solve_dc(); out[s]=cur()
        print(s*dv, r['iterations'], out[s], flush=True)
    G={k:(out[1][k]-out[-1][k])/(2*dv) for k in out[1]}
    print('dv',dv,'G',G, flush=True); res[dv]=G
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
# AC small signal
for f in (1.0,):
    ds.solve(type='ac', frequency=f)
    print('AC V1.I', f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
