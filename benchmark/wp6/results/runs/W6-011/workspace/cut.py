import numpy as np, devsim as ds
import mosfet_deck as M
M.build()
D,R=M.DEV,M.SI
x=np.array(ds.get_node_model_values(device=D,region=R,name='x'))
en=np.array(ds.get_element_node_list(device=D,region=R)) if False else None
e2n=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
n0=np.array(ds.get_edge_model_values(device=D,region=R,name='edge_index')) if False else None
ds.edge_from_node_model(device=D,region=R,node_model='x')
x0=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n0')); x1=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n1'))
def cuts():
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))
    out=[]
    for xc in (-0.2e-4, 0.5e-4, 1.0e-4, 1.5e-4, 2.2e-4):
        a=(x0<xc)&(x1>xc); b=(x1<xc)&(x0>xc)
        out.append(np.sum(J[a]*e2n[a])-np.sum(J[b]*e2n[b]))
    return np.array(out)
def I1(): return ds.get_circuit_node_value(solution='dcop', node='V1.I')
print('V0 cuts',cuts(), I1())
r={}
for v in (1e-4,-1e-4):
    ds.circuit_alter(name='V1',value=v); M.solve_dc(); r[v]=(cuts(),I1()); print(v, r[v], flush=True)
h=1e-4
print('G cuts', (r[h][0]-r[-h][0])/(2*h), 'G V1.I', (r[h][1]-r[-h][1])/(2*h))
