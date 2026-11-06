import json, numpy as np
import mosfet_deck as M
import devsim as ds
M.build()
D,R='mos','si'
ds.edge_from_node_model(device=D, region=R, node_model='x')
ds.edge_from_node_model(device=D, region=R, node_model='y')
x0=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n0')); x1=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n1'))
ec=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
xn=np.array(M.MESH['x_cm'])
cuts=[xn[i] for i in (20, 79, 150, 219, 300, 358, 420)]
def cutI():
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))
    out=[]
    for xc in cuts:
        # edges crossing x = xc + tiny (current from n0 to n1 direction)
        xm=xc+1e-12
        m=(np.minimum(x0,x1)<xm)&(np.maximum(x0,x1)>xm)
        s=np.sign(x1-x0)[m]   # orient to +x
        out.append(float(np.sum((J*ec)[m]*s)))
    return out
res={}
for dv in (1e-4,):
    o=[]
    for s in (1,-1):
        ds.circuit_alter(name='V1', value=s*dv); M.solve_dc()
        o.append((cutI(), ds.get_circuit_node_value(solution='dcop', node='V1.I')))
    G=[(a-b)/(2*dv) for a,b in zip(o[0][0],o[1][0])]
    res['cut_G']=G; res['V1_G']=(o[0][1]-o[1][1])/(2*dv)
    print(res, flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for f in (0.0, 1.0):
    ds.solve(type='ac', frequency=f, solver_type='direct')
    res[f'ac_{f}']=(ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
    print(res, flush=True)
json.dump(res, open('g_cut.json','w'), indent=1)
