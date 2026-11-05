import sys, time
import devsim as ds
import mosfet_deck as M
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    Id=ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation')
    Is=ds.get_contact_current(device=M.DEV, contact='source', equation='ElectronContinuityEquation')
    Ib=ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation')
    Iv=ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return Id,Is,Ib,Iv
print('V=0', cur(), flush=True)
res={}
for dv in [float(a) for a in sys.argv[1:]]:
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc(); res[dv]=cur(); print(dv, r['iterations'], res[dv], flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in [1e-3, 1.0]:
    ds.solve(type='ac', frequency=f)
    print('AC', f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'), flush=True)
