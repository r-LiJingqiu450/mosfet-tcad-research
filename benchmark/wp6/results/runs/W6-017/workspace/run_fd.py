import time, devsim as ds, mosfet_deck as M, numpy as np
t=time.time(); info=M.build(); print(info, time.time()-t)
def cur():
    i = {c: ds.get_contact_current(device='mos', contact=c, equation='ElectronContinuityEquation') for c in ('source','drain','body')}
    i['V1'] = ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return i
I0=cur(); print('V=0', I0)
res={}
for dv in (1e-3, -1e-3, 1e-4, -1e-4, 1e-5, -1e-5):
    ds.circuit_alter(name='V1', value=dv); M.solve_dc()
    res[dv]=cur(); print(dv, res[dv])
for dv in (1e-3,1e-4,1e-5):
    for k in ('drain','V1','source'):
        print(dv,k,(res[dv][k]-res[-dv][k])/(2*dv))
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in (1e-3, 1.0):
    ds.solve(type='ac', frequency=f)
    print('ac', f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'))
