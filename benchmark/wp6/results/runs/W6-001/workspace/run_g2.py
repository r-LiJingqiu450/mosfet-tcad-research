import numpy as np
import mosfet_deck as M
import devsim as ds
M.build()
D,R=M.DEV,M.SI
x,y=M.coords(R)
psi=np.array(ds.get_node_model_values(device=D,region=R,name='Potential'))
n=np.array(ds.get_node_model_values(device=D,region=R,name='Electrons'))
phin=psi-M.VT*np.log(n/M.NI)
print('phin range at 0 bias', phin.min(), phin.max())
en=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
n0=np.array(ds.get_edge_model_values(device=D,region=R,name='@n0')) if False else None
ds.edge_from_node_model(device=D,region=R,node_model='x'); ds.edge_from_node_model(device=D,region=R,node_model='y')
ds.edge_model(device=D,region=R,name='ex0',equation='x@n0'); ds.edge_model(device=D,region=R,name='ex1',equation='x@n1')
ds.edge_model(device=D,region=R,name='ey0',equation='y@n0'); ds.edge_model(device=D,region=R,name='ey1',equation='y@n1')
ex0,ex1,ey0,ey1=[np.array(ds.get_edge_model_values(device=D,region=R,name=k)) for k in ('ex0','ex1','ey0','ey1')]
xs=np.unique(x)
def cut_current(xc):
    # edges crossing x=xc
    m=((ex0<xc)&(ex1>xc))|((ex1<xc)&(ex0>xc))
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))
    sgn=np.sign(ex1-ex0)  # current flows n0->n1 positive? check
    return np.sum((J*en*sgn)[m])
cuts=[xs[100]+1e-9, 1e-4+1e-9, xs[300]+1e-9]
I0=[cut_current(c) for c in cuts]; print('cuts0',I0)
ds.solve(type='ac', frequency=1.0)
print('ac V1.I', ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'))
out=[]
for v in (1e-4,-1e-4):
    ds.circuit_alter(name='V1',value=v); M.solve_dc()
    out.append(([cut_current(c) for c in cuts], ds.get_circuit_node_value(node='V1.I', solution='dcop')))
    print(v,out[-1])
print('G cuts', [(a-b)/2e-4 for a,b in zip(out[0][0],out[1][0])], 'G V1', -(out[0][1]-out[1][1])/2e-4)
