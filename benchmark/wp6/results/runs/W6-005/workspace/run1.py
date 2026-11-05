import time, devsim as ds, mosfet_deck as M
t=time.time(); info=M.build(); print(info, time.time()-t)
def I(): 
    return (ds.get_contact_current(device='mos',contact='drain',equation='ElectronContinuityEquation'),
            ds.get_contact_current(device='mos',contact='source',equation='ElectronContinuityEquation'),
            ds.get_contact_current(device='mos',contact='body',equation='ElectronContinuityEquation'),
            ds.get_circuit_node_value(solution='dcop',node='V1.I'))
print('I0',I())
# small-signal AC
for f in (1e-3, 1.0):
    ds.solve(type='ac', frequency=f)
    print('ac',f, ds.get_circuit_node_value(solution='ssac_real',node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag',node='V1.I'))
res={}
for dv in (1e-4,-1e-4,1e-3,-1e-3):
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc(); res[dv]=I(); print(dv, r['iterations'], res[dv])
ds.circuit_alter(name='V1', value=0.0); M.solve_dc(); print('back I0', I())
for d in (1e-4,1e-3):
    print('central',d,[(res[d][k]-res[-d][k])/(2*d) for k in range(4)])
