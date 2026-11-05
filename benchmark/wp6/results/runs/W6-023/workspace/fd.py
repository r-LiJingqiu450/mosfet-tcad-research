import time, devsim as ds, mosfet_deck as M
t=time.time(); print(M.build(), time.time()-t, flush=True)
def cur():
    return (ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation'),
            ds.get_circuit_node_value(node='V1.I'),
            ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation'),
            ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation'))
print('V=0', cur(), flush=True)
res={}
for dv in (1e-3, 2e-3, 4e-4):
    for s in (1,-1):
        ds.circuit_alter(name='V1', value=s*dv); r=M.solve_dc()
        res[s*dv]=cur(); print(s*dv, r['iterations'], res[s*dv], flush=True)
    for k in range(2):
        print('dv',dv,'G[%d]'%k, (res[dv][k]-res[-dv][k])/(2*dv), flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in (0.0, 1.0):
    ds.solve(type='ac', frequency=f)
    print('AC f',f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'), ds.get_circuit_solution_list())
