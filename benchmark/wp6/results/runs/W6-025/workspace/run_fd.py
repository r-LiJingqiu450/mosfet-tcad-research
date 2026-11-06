import time, json
import mosfet_deck as M
import devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    ie=ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation')
    isrc=ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation')
    ib=ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation')
    ic=ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return ie, isrc, ib, ic
res={}
res[0.0]=cur(); print(0.0, res[0.0], flush=True)
for v in [1e-4,-1e-4,1e-3,-1e-3]:
    ds.circuit_alter(name='V1', value=v); r=M.solve_dc()
    res[v]=cur(); print(v, r['iterations'], res[v], flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in [1.0, 1e3]:
    ds.solve(type='ac', frequency=f)
    print('ac', f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'), flush=True)
for a,b in [(1e-4,-1e-4),(1e-3,-1e-3)]:
    for k in range(4):
        print(a, k, (res[a][k]-res[b][k])/(a-b))
