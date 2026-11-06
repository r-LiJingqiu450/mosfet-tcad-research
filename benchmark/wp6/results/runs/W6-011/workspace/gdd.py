import sys, time
import devsim as ds
import mosfet_deck as M
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    Ic = ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation')
    Is = ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation')
    Ib = ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation')
    Iv = ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return Ic, Is, Ib, Iv
I0=cur(); print('V=0', I0, flush=True)
res={}
for h in [1e-3, 1e-4, 1e-5]:
    r=[]
    for v in (h,-h):
        ds.circuit_alter(name='V1', value=v); M.solve_dc(); r.append(cur()); print(v, r[-1], flush=True)
    g=[(a-b)/(2*h) for a,b in zip(*r)]
    print('h',h,'G (drain, source, body, V1.I):', g, flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in (1e-3, 1.0):
    ds.solve(type='ac', frequency=f)
    print('AC f',f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'), flush=True)
