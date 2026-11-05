import time, json
import mosfet_deck as M
import devsim as ds
t=time.time()
print(M.build(), time.time()-t)
def cur():
    d={c: ds.get_contact_current(device='mos', contact=c, equation='ElectronContinuityEquation') for c in ('source','drain','body')}
    d['V1']=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return d
I0=cur(); print('V=0', I0)
res={}
for dv in (1e-3, 1e-4, 1e-5):
    out=[]
    for s in (+1,-1):
        ds.circuit_alter(name='V1', value=s*dv); r=M.solve_dc(); out.append(cur()); print(s*dv, r['iterations'], out[-1], flush=True)
    G={k:(out[0][k]-out[1][k])/(2*dv) for k in out[0]}
    print('dv',dv,'G',G, flush=True); res[dv]=G
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
json.dump({str(k):v for k,v in res.items()}, open('g_fd.json','w'), indent=1)
