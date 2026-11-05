import numpy as np, devsim as ds, mosfet_deck as M
M.build()
D,R=M.DEV,M.SI
ds.edge_from_node_model(device=D, region=R, node_model='x')
def e(n): return np.array(ds.get_edge_model_values(device=D, region=R, name=n))
x0,x1=e('x@n0'),e('x@n1')
xs=np.array(M.MESH['x_cm'])
cuts=[(xs[i]+xs[i+1])/2 for i in (5, 60, 150, 219, 300, 380, 430)]
def fluxes():
    J=e('ElectronCurrent')*e('EdgeCouple')
    out=[]
    for c in cuts:
        m=(np.minimum(x0,x1)<c)&(np.maximum(x0,x1)>c)
        s=np.sign(x1-x0)[m]
        out.append(np.sum(J[m]*s))
    return np.array(out)
F0=fluxes(); print('cuts', cuts); print('F0', F0)
res={}
for d in (1e-6,-1e-6,1e-4,-1e-4):
    ds.circuit_alter(name='V1', value=d); M.solve_dc()
    F=fluxes(); I=ds.get_circuit_node_value(solution='dcop', node='V1.I')
    res[d]=(F,I)
    print('d',d,'(F-F0)/d', (F-F0)/d, '-V1.I/d', -I/d, flush=True)
for h in (1e-6,1e-4):
    print('central h',h,'cuts',(res[h][0]-res[-h][0])/(2*h),'circ',-(res[h][1]-res[-h][1])/(2*h))
