import numpy as np, mosfet_deck as M, devsim as ds
D,S=M.DEV,M.SI
M.build()
x,y=M.coords(S)
def ev(name): return np.array(ds.get_edge_model_values(device=D,region=S,name=name))
n0=np.array(ds.get_edge_model_values(device=D,region=S,name='@n0') if False else [])
ds.edge_from_node_model(device=D,region=S,node_model='x'); ds.edge_from_node_model(device=D,region=S,node_model='y')
def cuts(xc_list):
    J=ev('ElectronCurrent'); C=ev('EdgeCouple'); x0=ev('x@n0'); x1=ev('x@n1')
    out=[]
    for xc in xc_list:
        m1=(x0<xc)&(x1>xc); m2=(x1<xc)&(x0>xc)
        out.append(np.sum(J[m1]*C[m1])-np.sum(J[m2]*C[m2]))
    return np.array(out)
xcs=np.array([-0.2,-0.05,0.1,0.3,0.5,0.7,0.9,1.05,1.2])*1e-4
xcs=xcs+1e-9
res={}
for v in [0,1e-4,-1e-4]:
    ds.circuit_alter(name='V1',value=v); M.solve_dc()
    res[v]=cuts(xcs); print(v, res[v], ds.get_circuit_node_value(node='V1.I',solution='dcop'))
print('G cuts', (res[1e-4]-res[-1e-4])/2e-4)
ds.circuit_alter(name='V1',value=0); M.solve_dc()
psi=np.array(ds.get_node_model_values(device=D,region=S,name='Potential'))
n=np.array(ds.get_node_model_values(device=D,region=S,name='Electrons'))
xu=np.unique(x)
for xv in [-0.1e-4,0,0.05e-4,0.25e-4,0.5e-4,0.75e-4,1e-4]:
    i=np.argmin(abs(xu-xv)); m=(x==xu[i])&(y<1e-12)&(y>-1e-12)
    print('surface x=%.4g psi=%.4f n=%.3e'%(xu[i]*1e4, psi[m][0], n[m][0]))
