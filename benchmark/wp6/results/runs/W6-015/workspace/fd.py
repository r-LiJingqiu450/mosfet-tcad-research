import time, json
import devsim as ds
import mosfet_deck as M
t=time.time()
print(M.build(), time.time()-t)
def cur():
    Id = ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation')
    Is = ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation')
    Ib = ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation')
    Iv = ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return Id, Is, Ib, Iv
print('V=0', cur())
# AC small signal
ds.solve(type='ac', frequency=1.0)
print('AC V1.I', ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'))
res={}
for dv in (1e-3, 1e-4, 1e-5):
    out=[]
    for s in (dv,-dv):
        ds.circuit_alter(name='V1', value=s); M.solve_dc(); out.append(cur())
    ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
    g=[(a-b)/(2*dv) for a,b in zip(out[0],out[1])]
    print(dv, out, 'G', g, flush=True)
