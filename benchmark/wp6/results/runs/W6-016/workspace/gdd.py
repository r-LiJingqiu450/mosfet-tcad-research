import time, sys
import mosfet_deck as M
import devsim as ds
t=time.time()
info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    Ic=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    Id=ds.get_contact_current(device='mos', contact='drain', equation='ElectronContinuityEquation')
    Is=ds.get_contact_current(device='mos', contact='source', equation='ElectronContinuityEquation')
    Ib=ds.get_contact_current(device='mos', contact='body', equation='ElectronContinuityEquation')
    return Ic,Id,Is,Ib
print('V=0', cur(), flush=True)
# AC small-signal at low frequency
for f in (1.0, 1e3):
    ds.solve(type='ac', frequency=f)
    print('AC f',f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'), flush=True)
res={}
for dv in (1e-3,1e-4,1e-5):
    for s in (+1,-1):
        ds.circuit_alter(name='V1', value=s*dv); M.solve_dc(); res[s*dv]=cur(); print(s*dv,res[s*dv],flush=True)
    print('dv',dv,'G central (circuit)',(res[dv][0]-res[-dv][0])/(2*dv),'(drain contact)',(res[dv][1]-res[-dv][1])/(2*dv),flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc(); print('back to 0', cur())
