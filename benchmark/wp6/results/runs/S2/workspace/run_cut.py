import numpy as np
import mosfet_deck as M
import devsim as ds
M.build()
D,R='mos','si'
def edges():
    n0=np.array(ds.get_edge_model_values(device=D,region=R,name='@n0') if False else [],int)
x=np.array(ds.get_node_model_values(device=D,region=R,name='x'))
y=np.array(ds.get_node_model_values(device=D,region=R,name='y'))
el=np.array(ds.get_element_node_list(device=D,region=R)) if False else None
ec=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
# edge node indices
ds.edge_from_node_model(device=D,region=R,node_model='x'); ds.edge_from_node_model(device=D,region=R,node_model='y')
x0=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n0')); x1=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n1'))
y0=np.array(ds.get_edge_model_values(device=D,region=R,name='y@n0')); y1=np.array(ds.get_edge_model_values(device=D,region=R,name='y@n1'))
xs=np.unique(x)
def cutI():
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))
    out=[]
    for xc in [xs[10],xs[60],xs[130],xs[200],xs[250]]:
        xm=xc+1e-12
        cross=((x0<xm)&(x1>xm))|((x1<xm)&(x0>xm))
        s=np.where(x1>x0,1,-1)  # J positive n0->n1 ? sign convention: current from n0 to n1
        out.append(np.sum((J*ec*s)[cross]))
    return np.array(out)
I0=cutI(); print('I0',I0)
for dv in (1e-3,-1e-3,1e-4):
    ds.circuit_alter(name='V1',value=dv); M.solve_dc()
    I=cutI(); print('RES',dv,'cut I',I,'G',(I-I0)/dv,'V1.I/dv',ds.get_circuit_node_value(node='V1.I',solution='dcop')/dv,flush=True)
n=np.array(ds.get_node_model_values(device=D,region=R,name='Electrons'))
m=(np.abs(x-xs[130])<1e-12); print('n at mid col (top):',n[m][np.argsort(y[m])][:5])
