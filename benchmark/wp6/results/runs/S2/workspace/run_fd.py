import time, json
import mosfet_deck as M
import devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    return (ds.get_contact_current(device='mos',contact='drain',equation='ElectronContinuityEquation'),
            ds.get_contact_current(device='mos',contact='source',equation='ElectronContinuityEquation'),
            ds.get_contact_current(device='mos',contact='body',equation='ElectronContinuityEquation'),
            ds.get_circuit_node_value(node='V1.I', solution='dcop'))
print('V=0', cur(), flush=True)
res={}
for dv in (1e-3,-1e-3,1e-4,-1e-4,1e-5,-1e-5):
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc(); c=cur(); res[dv]=c
    print(dv, r['iterations'], c, c[0]/dv, c[3]/dv, flush=True)
for d in (1e-3,1e-4,1e-5):
    print('central', d, (res[d][0]-res[-d][0])/(2*d), (res[d][3]-res[-d][3])/(2*d))
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in (1e-3, 1.0):
    ds.solve(type='ac', frequency=f)
    print('ac', f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'))
