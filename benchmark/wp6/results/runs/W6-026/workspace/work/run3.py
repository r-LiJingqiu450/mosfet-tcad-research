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
# one Newton step test
h=1e-5
ds.circuit_alter(name='V1', value=h); M._solve_one(); print('ONESTEP', cur()[2]/h, flush=True)
M.solve_dc(); ip=cur(); print('R',h,ip, time.time()-t, flush=True)
ds.circuit_alter(name='V1', value=-h); M.solve_dc(); im=cur(); print('R',-h,im, flush=True)
print('FD', [(a-b)/(2*h) for a,b in zip(ip,im)], flush=True)
