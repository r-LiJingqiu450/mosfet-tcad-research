import sys, time
sys.path.insert(0, '/home/claude/wp6_runs/W6-031')
import numpy as np, devsim as ds
import mosfet_deck as M
EXT = sys.argv[1] == 'ext'
if EXT:
    for p in ('extended_solver','extended_model','extended_equation'):
        ds.set_parameter(name=p, value=True)
t=time.time(); info=M.build(); print(info, time.time()-t, flush=True)
D,R=M.DEV,M.SI
ds.edge_from_node_model(device=D, region=R, node_model='x')
ds.edge_from_node_model(device=D, region=R, node_model='y')
def cut(xc):
    x0=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n0'))
    x1=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n1'))
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))
    C=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
    m1=(x0<xc)&(x1>xc); m2=(x1<xc)&(x0>xc)
    # ElectronCurrent is the flux from n0 to n1 direction (DEVSIM convention); return signed sum
    return float(np.sum(J[m1]*C[m1])-np.sum(J[m2]*C[m2]))
def cur():
    out={c:ds.get_contact_current(device=D, contact=c, equation='ElectronContinuityEquation') for c in ('source','drain','body')}
    out['V1.I']=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    for xc in (0.25e-4,0.5e-4,0.75e-4): out[f'cut{xc}']=cut(xc)
    return out
print('V=0', cur(), flush=True)
for f in (1e-3,1.0):
    ds.solve(type='ac', frequency=f)
    print('AC', f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'), flush=True)
res={}
for v in (1e-4,-1e-4,1e-3,-1e-3):
    ds.circuit_alter(name='V1', value=v); h=M.solve_dc()
    res[v]=cur(); print(v, h['iterations'], res[v], flush=True)
for h in (1e-4,1e-3):
    print(h, {k:(res[h][k]-res[-h][k])/(2*h) for k in res[h]})
