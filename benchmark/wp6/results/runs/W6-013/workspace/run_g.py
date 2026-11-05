import sys, time, numpy as np
import mosfet_deck as M, devsim as ds
t=time.time(); print(M.build(), time.time()-t, flush=True)
def cur():
    Id = ds.get_contact_current(device='mos', contact='drain', equation='ElectronContinuityEquation')
    Is = ds.get_contact_current(device='mos', contact='source', equation='ElectronContinuityEquation')
    Ib = ds.get_contact_current(device='mos', contact='body', equation='ElectronContinuityEquation')
    try: Iv = ds.get_circuit_node_value(solution='dcop', node='V1.I')
    except Exception as e: Iv = float('nan')
    return Id, Is, Ib, Iv
print('V=0', cur(), flush=True)
# AC small-signal
ds.solve(type='ac', frequency=1e-3)
try:
    print('AC V1.I', ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'), flush=True)
except Exception as e: print('ac err', e)
res={}
for dv in [float(a) for a in sys.argv[1:]]:
    for s in (+1,-1):
        ds.circuit_alter(name='V1', value=s*dv); r=M.solve_dc(); res[(dv,s)]=cur()
        print(s*dv, r['iterations'], res[(dv,s)], flush=True)
    G=[(res[(dv,1)][k]-res[(dv,-1)][k])/(2*dv) for k in range(4)]
    print('dv',dv,'G (Id,Is,Ib,IV1)',G, flush=True)
ds.circuit_alter(name='V1', value=0.); M.solve_dc()
