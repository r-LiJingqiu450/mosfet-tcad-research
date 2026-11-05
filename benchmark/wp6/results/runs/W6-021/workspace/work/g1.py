import sys, time
sys.path.insert(0, '.')
import mosfet_deck as M
import devsim as ds
t=time.time()
info = M.build(); print(info, time.time()-t)
def cur():
    I = ds.get_circuit_node_value(solution='dcop', node='V1.I')
    Ie = ds.get_contact_current(device='mos', contact='drain', equation='ElectronContinuityEquation')
    Is = ds.get_contact_current(device='mos', contact='source', equation='ElectronContinuityEquation')
    Ib = ds.get_contact_current(device='mos', contact='body', equation='ElectronContinuityEquation')
    return I, Ie, Is, Ib
print('V=0', cur())
res={}
for dv in [1e-3, -1e-3, 1e-4, -1e-4, 1e-5,-1e-5]:
    ds.circuit_alter(name='V1', value=dv); M.solve_dc()
    res[dv]=cur(); print(dv, res[dv], flush=True)
for d in [1e-3,1e-4,1e-5]:
    print(d, [(res[d][k]-res[-d][k])/(2*d) for k in range(4)])
ds.circuit_alter(name='V1', value=0.); M.solve_dc()
for f in [0.0, 1.0, 1e3]:
    ds.solve(type='ac', frequency=f)
    print('ac', f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
