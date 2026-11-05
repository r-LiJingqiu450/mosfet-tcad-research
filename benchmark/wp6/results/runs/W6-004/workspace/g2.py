import sys, numpy as np, devsim as ds, mosfet_deck as M
D,S=M.DEV,M.SI
if len(sys.argv)>1 and sys.argv[1]=='ext':
    for p in ('extended_solver','extended_model','extended_equation'): ds.set_parameter(name=p, value=True)
print(M.build(), flush=True)
for nm in ('x','y'):
    ds.edge_from_node_model(device=D, region=S, node_model=nm)
def ev(n): return np.array(ds.get_edge_model_values(device=D, region=S, name=n))
x0,x1=ev('x@n0'),ev('x@n1')
cuts=[0.25e-4,0.5e-4,1e-4,1.5e-4,1.75e-4]
def meas():
    J=ev('ElectronCurrent'); C=ev('EdgeCouple')
    out=[]
    for xc in cuts:
        a=(x0<xc)&(x1>xc); b=(x1<xc)&(x0>xc)
        out.append(np.sum(J[a]*C[a])-np.sum(J[b]*C[b]))
    return dict(V1I=ds.get_circuit_node_value(node='V1.I', solution='dcop'),
                Id=ds.get_contact_current(device=D, contact='drain', equation='ElectronContinuityEquation'),
                Ib=ds.get_contact_current(device=D, contact='body', equation='ElectronContinuityEquation'),
                cut=np.array(out))
r={0:meas()}; print(0, r[0])
for v in [1e-4,-1e-4,1e-3,-1e-3]:
    ds.circuit_alter(name='V1', value=v); M.solve_dc(); r[v]=meas(); print(v, r[v], flush=True)
for d in [1e-4,1e-3]:
    print('G', d, {k:(r[d][k]-r[-d][k])/(2*d) for k in r[d]})
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in [1.0, 1e3]:
    ds.solve(type='ac', frequency=f)
    print('AC', f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'))
