import time, json
import devsim as ds
import mosfet_deck as M
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    r={}
    for c in ('drain','source','body'):
        r[c]=ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation')
    r['V1.I']=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    r['gate']=ds.get_contact_charge(device=M.DEV, contact='gate', equation='PotentialEquation')
    return r
res={0.0:cur()}; print(0.0,res[0.0],flush=True)
for v in (1e-3,-1e-3,2.5e-4,-2.5e-4):
    ds.circuit_alter(name='V1', value=v)
    h=M.solve_dc()
    res[v]=cur(); print(v,h['iterations'],res[v],flush=True)
for d in (1e-3,2.5e-4):
    for k in ('drain','source','body','V1.I'):
        print(d,k,'central G',(res[d][k]-res[-d][k])/(2*d),'fwd',(res[d][k]-res[0.0][k])/d)
json.dump({str(k):v for k,v in res.items()},open('fd_results.json','w'),indent=1)
