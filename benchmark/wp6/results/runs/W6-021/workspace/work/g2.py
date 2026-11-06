import sys
sys.path.insert(0, '.')
import numpy as np
import mosfet_deck as M
import devsim as ds
M.build()
D,R='mos','si'
ds.edge_from_node_model(device=D, region=R, node_model='node_index')
ds.edge_from_node_model(device=D, region=R, node_model='x')
ds.edge_from_node_model(device=D, region=R, node_model='y')
x,y=M.coords()
ymax=y.max()
i0=np.array(ds.get_edge_model_values(device=D,region=R,name='node_index@n0')).astype(int)
i1=np.array(ds.get_edge_model_values(device=D,region=R,name='node_index@n1')).astype(int)
ec=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
body = np.abs(y-ymax)<1e-15
def flux(xc):
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))*ec
    S=(x>xc)&~body
    a=S[i0]&~S[i1]; b=S[i1]&~S[i0]
    return J[a].sum()-J[b].sum()
out={}
for v in [0., 1e-3,-1e-3,1e-4,-1e-4]:
    ds.circuit_alter(name='V1', value=v); M.solve_dc()
    I=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    out[v]=[I]+[flux(xc) for xc in (0.2e-4,0.5e-4,0.8e-4)]
    print(v,out[v],flush=True)
for d in [1e-3,1e-4]:
    print(d, (np.array(out[d])-np.array(out[-d]))/(2*d))
