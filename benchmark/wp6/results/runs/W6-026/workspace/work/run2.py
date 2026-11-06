import sys, time; sys.path.insert(0,'.')
import devsim as ds, mosfet_deck as M
ext = len(sys.argv)>1 and sys.argv[1]=='ext'
if ext:
    for p in ('extended_solver','extended_model','extended_equation'):
        ds.set_parameter(name=p, value=True)
t=time.time(); info=M.build(); print('BUILD',info, time.time()-t, flush=True)
def cur():
    return (ds.get_contact_current(device='mos',contact='drain',equation='ElectronContinuityEquation'),
            ds.get_contact_current(device='mos',contact='source',equation='ElectronContinuityEquation'),
            ds.get_circuit_node_value(solution='dcop', node='V1.I'))
def ac():
    ds.solve(type='ac', frequency=1e-3)
    return ds.get_circuit_node_value(solution='ssac_real', node='V1.I')
print('R0', cur(), 'AC', ac(), flush=True)
for k in range(3):
    M._solve_one(); print('R0 extra', cur(), flush=True)
res={}
for h in (1e-5,1e-4,1e-3,3e-3):
    for s in (1,-1):
        ds.circuit_alter(name='V1', value=s*h); M.solve_dc(); res[s*h]=cur()
        print('R', s*h, res[s*h], 'AC', ac(), flush=True)
    for i,nm in enumerate(('drain','source','V1')):
        print('FD', h, nm, (res[h][i]-res[-h][i])/(2*h), flush=True)
