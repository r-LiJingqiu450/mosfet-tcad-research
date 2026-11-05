import numpy as np, devsim as ds, mosfet_deck as M
M.build()
x,y=M.coords()
en=np.array(ds.get_edge_model_values(device='mos',region='si',name='EdgeNodes') if False else [])
e0=np.array(ds.get_edge_model_values(device='mos',region='si',name='edge_index') if False else [])
# edge node indices
ds.edge_from_node_model(device='mos',region='si',node_model='x')
x0=np.array(ds.get_edge_model_values(device='mos',region='si',name='x@n0'))
x1=np.array(ds.get_edge_model_values(device='mos',region='si',name='x@n1'))
cpl=np.array(ds.get_edge_model_values(device='mos',region='si',name='EdgeCouple'))
def cut(xc):
    J=np.array(ds.get_edge_model_values(device='mos',region='si',name='ElectronCurrent'))
    m1=(x0<xc)&(x1>xc); m2=(x1<xc)&(x0>xc)
    return np.sum(J[m1]*cpl[m1])-np.sum(J[m2]*cpl[m2])
xcs=[0.5e-4,1.0e-4,1.5e-4]
r={}
for dv in [1e-3,-1e-3]:
    ds.circuit_alter(name='V1', value=dv); M.solve_dc()
    r[dv]=[ds.get_circuit_node_value(solution='dcop', node='V1.I')]+[cut(xc) for xc in xcs]
print('CUT G', [(a-b)/2e-3 for a,b in zip(r[1e-3],r[-1e-3])])
