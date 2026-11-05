import numpy as np, devsim as ds, mosfet_deck as M
M.build()
D,R=M.DEV,M.SI
x=np.array(ds.get_node_model_values(device=D,region=R,name='x'))
xs=np.unique(x)
e0=np.array(ds.get_edge_model_values(device=D,region=R,name='node0_index' if False else 'EdgeCouple'))
# edge node indices
ei=np.array(ds.get_element_node_list(device=D,region=R)) if False else None
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
en=np.stack([np.array(ds.get_edge_model_values(device=D,region=R,name='node_index@n%d'%k)).astype(int) for k in (0,1)],1)
print('edge list', None if en is None else en.shape)
def cuts():
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))
    out=[]
    for xc in (xs[40], xs[150], xs[184], xs[250], xs[330]):
        a=x[en[:,0]]<=xc; b=x[en[:,1]]<=xc
        cross=a!=b
        s=np.where(a[cross],1.0,-1.0)  # current from n0 to n1 when n0 left
        out.append(np.sum(s*J[cross]*e0[cross]))
    return np.array(out)
I={}
for v in (1e-4,-1e-4):
    ds.circuit_alter(name='V1',value=v); M.solve_dc()
    I[v]=cuts(); print(v, I[v], ds.get_circuit_node_value(node='V1.I',solution='dcop'))
print('G cuts', (I[1e-4]-I[-1e-4])/2e-4)
