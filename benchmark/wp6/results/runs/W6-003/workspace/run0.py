import time, devsim as ds, mosfet_deck as M, numpy as np
t=time.time(); info=M.build(); print(info, time.time()-t)
def cur():
    Id=ds.get_contact_current(device='mos',contact='drain',equation='ElectronContinuityEquation')
    Is=ds.get_contact_current(device='mos',contact='source',equation='ElectronContinuityEquation')
    Ib=ds.get_contact_current(device='mos',contact='body',equation='ElectronContinuityEquation')
    Ic=ds.get_circuit_node_value(solution='dcop',node='V1.I')
    return Id,Is,Ib,Ic
print('V=0',cur())
for dv in [1e-3,-1e-3]:
    ds.circuit_alter(name='V1',value=dv); r=M.solve_dc(); print(dv,r['iterations'],cur())
