import numpy as np, devsim as ds, mosfet_deck as M
M.build()
D,R=M.DEV,M.SI
for v in ('x','y','Electrons'):
    ds.edge_from_node_model(device=D,region=R,node_model=v)
xs=np.array(M.MESH['x_cm'])
cuts=[0.5*(xs[i]+xs[i+1]) for i in (40,100,183,250,330)]
def flux():
    ds.edge_model(device=D,region=R,name='Jtmp',equation='ElectronCurrent*EdgeCouple')
    g=lambda n: np.array(ds.get_edge_model_values(device=D,region=R,name=n))
    F,x0,x1=g('Jtmp'),g('x@n0'),g('x@n1')
    out=[]
    for xc in cuts:
        m=(np.minimum(x0,x1)<xc)&(np.maximum(x0,x1)>xc)
        # ElectronCurrent is flux from n0 to n1? sign: orient to +x
        s=np.where(x1>x0,1.0,-1.0)
        out.append(np.sum(F[m]*s[m]))
    Ic=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return np.array(out+[Ic])
f0=flux(); print('V=0',f0)
res={}
for dv in (3e-4,1e-4,-1e-4,-3e-4):
    ds.circuit_alter(name='V1',value=dv); M.solve_dc(); res[dv]=flux(); print(dv,res[dv],flush=True)
for dv in (3e-4,1e-4):
    print('central',dv,(res[dv]-res[-dv])/(2*dv))
    print('fwd',dv,(res[dv]-f0)/dv)
