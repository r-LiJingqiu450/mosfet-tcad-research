import time, devsim as ds, mosfet_deck as M
t=time.time(); info=M.build(); print(info, time.time()-t, flush=True)
def cur():
    d={c: ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation') for c in ('source','drain','body')}
    d['V1.I']=ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return d
I0=cur(); print('V=0', I0, flush=True)
ds.solve(type='ac', frequency=1.0)
print('AC f=1 V1.I real/imag', ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'), flush=True)
res={}
for dv in (1e-3, -1e-3, 1e-4, -1e-4, 1e-5, -1e-5):
    ds.circuit_alter(name='V1', value=dv); M.solve_dc(); res[dv]=cur(); print(dv, res[dv], flush=True)
for dv in (1e-3,1e-4,1e-5):
    for k in ('drain','V1.I','source','body'):
        print(dv, k, 'central G', (res[dv][k]-res[-dv][k])/(2*dv))
