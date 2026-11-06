import time, numpy as np, devsim as ds, mosfet_deck as M
t=time.time(); info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    return (ds.get_circuit_node_value(solution='dcop', node='V1.I'),
            ds.get_contact_current(device='mos', contact='drain', equation='ElectronContinuityEquation'),
            ds.get_contact_current(device='mos', contact='source', equation='ElectronContinuityEquation'),
            ds.get_contact_current(device='mos', contact='body', equation='ElectronContinuityEquation'))
print('V=0', cur(), flush=True)
res={}
for dv in [1e-3,-1e-3,1e-4,-1e-4,1e-5,-1e-5]:
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc(); res[dv]=cur(); print(dv, r['iterations'], res[dv], flush=True)
for dv in [1e-3,1e-4,1e-5]:
    print('G central', dv, [(a-b)/(2*dv) for a,b in zip(res[dv],res[-dv])])
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in [1e-3, 1.0, 1e3]:
    ds.solve(type='ac', frequency=f)
    print('AC', f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'), flush=True)
