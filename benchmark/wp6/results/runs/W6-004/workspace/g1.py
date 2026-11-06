import time, devsim as ds, mosfet_deck as M
t=time.time()
print(M.build(), time.time()-t)
def cur():
    i_e = ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation')
    i_c = ds.get_circuit_node_value(node='V1.I', solution='dcop')
    i_s = ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation')
    i_b = ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation')
    return i_e, i_c, i_s, i_b
print('V=0', cur())
res={}
for v in [1e-4,-1e-4,1e-3,-1e-3, 1e-5,-1e-5]:
    ds.circuit_alter(name='V1', value=v); r=M.solve_dc(); c=cur(); res[v]=c
    print(v, r['iterations'], c, flush=True)
for d in [1e-5,1e-4,1e-3]:
    print(d, [(res[d][k]-res[-d][k])/(2*d) for k in range(4)])
