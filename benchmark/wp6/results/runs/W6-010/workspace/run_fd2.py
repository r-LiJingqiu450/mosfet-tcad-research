import sys, json
import devsim as ds
import mosfet_deck as M
if len(sys.argv)>1 and sys.argv[1]=='ext':
    for p in ('extended_solver','extended_model','extended_equation'): ds.set_parameter(name=p, value=True)
M.build()
def cur():
    return [ds.get_circuit_node_value(solution='dcop', node='V1.I'),
            ds.get_contact_current(device='mos',contact='drain',equation='ElectronContinuityEquation'),
            ds.get_contact_current(device='mos',contact='source',equation='ElectronContinuityEquation')]
def ac():
    ds.solve(type='ac', frequency=0.0)
    return ds.get_circuit_node_value(solution='ssac_real', node='V1.I')
out={}
out['0']=cur()+[ac()]
print('RES 0', out['0'], flush=True)
for dv in [1e-5,-1e-5,1e-4,-1e-4,1e-3,-1e-3,3e-3,-3e-3,1e-2,-1e-2]:
    ds.circuit_alter(name='V1', value=dv); M.solve_dc()
    out[str(dv)]=cur()+[ac()]
    print('RES',dv,out[str(dv)],flush=True)
json.dump(out,open(f'fd2_{sys.argv[1] if len(sys.argv)>1 else "dbl"}.json','w'),indent=1)
