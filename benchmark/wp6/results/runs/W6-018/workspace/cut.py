import numpy as np, devsim as ds, mosfet_deck as M, sys
info=M.build()
for q in ('x','y'):
    ds.edge_from_node_model(device=M.DEV, region=M.SI, node_model=q)
g=lambda n: np.array(ds.get_edge_model_values(device=M.DEV, region=M.SI, name=n))
x0,x1,y0,y1=g('x@n0'),g('x@n1'),g('y@n0'),g('y@n1')
EC=g('EdgeCouple')
xs=np.array(M.MESH['x_cm'])
cuts=[(xs[i]+xs[i+1])/2 for i in (40,100,219,300,400)]
def cutI():
    J=g('ElectronCurrent'); out=[]
    for c in cuts:
        s=np.where((x0<c)&(x1>c),1.0,np.where((x1<c)&(x0>c),-1.0,0.0))
        out.append(np.sum(s*J*EC))  # flux from n0 to n1 direction sign as defined
    # recompute with longdouble from psi, n
    return np.array(out)
def recompute(c):
    psi0,psi1,n0,n1=g('Potential@n0'),g('Potential@n1'),g('Electrons@n0'),g('Electrons@n1')
    L=g('EdgeLength')
    m=((x0<c)&(x1>c))|((x1<c)&(x0>c))
    s=np.where(x0<c,1,-1)[m]
    ld=np.longdouble
    v=(psi0[m].astype(ld)-psi1[m].astype(ld))/ld(M.VT)
    B=lambda z: np.where(np.abs(z)<1e-8, 1-z/2, z/np.expm1(z))
    j=ld(M.Q*M.MU_N*M.VT)/L[m].astype(ld)*(n1[m].astype(ld)*B(v)-n0[m].astype(ld)*B(-v))
    return float(np.sum(s*j*EC[m].astype(ld)))
def all_I():
    return dict(V1=ds.get_circuit_node_value(node='V1.I', solution='dcop'),
                drain=ds.get_contact_current(device=M.DEV, contact='drain', equation='ElectronContinuityEquation'),
                cuts=cutI(), rec=recompute(cuts[2]))
R={0:all_I()}
print(0,R[0],flush=True)
for dv in (1e-3,-1e-3,2e-4,-2e-4):
    ds.circuit_alter(name='V1', value=dv); M.solve_dc(); R[dv]=all_I(); print(dv,R[dv],flush=True)
for dv in (1e-3,2e-4):
    print('dv',dv,'V1 G',(R[dv]['V1']-R[-dv]['V1'])/2/dv,'cut G',(R[dv]['cuts']-R[-dv]['cuts'])/2/dv,'rec G',(R[dv]['rec']-R[-dv]['rec'])/2/dv)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in (1e-3,1.0,1e3):
    ds.solve(type='ac', frequency=f)
    print('AC',f,ds.get_circuit_node_value(node='V1.I', solution='ssac_real'),ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'))
