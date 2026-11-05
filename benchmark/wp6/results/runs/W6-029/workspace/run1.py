import mosfet_deck as M, devsim as ds, time
t=time.time()
info=M.build(); print(info, time.time()-t)
def cur():
    Ic = ds.get_contact_current(device='mos', contact='drain', equation='ElectronContinuityEquation')
    Is = ds.get_contact_current(device='mos', contact='source', equation='ElectronContinuityEquation')
    Ib = ds.get_contact_current(device='mos', contact='body', equation='ElectronContinuityEquation')
    try: Iv = ds.get_circuit_node_value(solution='dcop', node='V1.I')
    except Exception as e: Iv=str(e)
    return Ic,Is,Ib,Iv
print('V=0', cur())
for v in [1e-3,-1e-3,1e-4,-1e-4]:
    ds.circuit_alter(name='V1', value=v); r=M.solve_dc(); c=cur(); print(v, r['iterations'], c, c[0]/v)
