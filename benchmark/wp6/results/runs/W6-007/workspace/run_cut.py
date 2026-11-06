import numpy as np, devsim as ds, mosfet_deck as M
M.build()
ds.edge_from_node_model(device=M.DEV, region=M.SI, node_model='x')
ds.edge_from_node_model(device=M.DEV, region=M.SI, node_model='y')
g=lambda n: np.array(ds.get_edge_model_values(device=M.DEV, region=M.SI, name=n))
x0,x1,y0,y1,ec=g('x@n0'),g('x@n1'),g('y@n0'),g('y@n1'),g('EdgeCouple')
xs=np.array(M.MESH['x_cm'])
icuts=[5,50,100,150,219,300,350,400,432]
def cuts():
    J=g('ElectronCurrent')*ec
    out=[]
    for i in icuts:
        xc=0.5*(xs[i]+xs[i+1])
        m=(np.minimum(x0,x1)<xc)&(np.maximum(x0,x1)>xc)
        s=np.sign(x1-x0)[m]   # current n0->n1 positive? orient along +x
        out.append(np.sum(J[m]*s))
    Ib=ds.get_contact_current(device=M.DEV, contact='body', equation='ElectronContinuityEquation')
    return np.array(out), Ib, ds.get_circuit_node_value(node='V1.I', solution='dcop')
np.set_printoptions(linewidth=250,precision=6)
base=cuts(); print('V=0', base)
r={}
for dv in [1e-4,-1e-4,1e-3,-1e-3]:
    ds.circuit_alter(name='V1', value=dv); M.solve_dc(); r[dv]=cuts()
for d in [1e-4,1e-3]:
    a,b=r[d],r[-d]
    print('dV',d,'G_cut',(a[0]-b[0])/(2*d),'G_body',(a[1]-b[1])/(2*d),'G_V1',(a[2]-b[2])/(2*d))
print('x cuts (um)', [0.5e4*(xs[i]+xs[i+1]) for i in icuts])
