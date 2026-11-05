import numpy as np, devsim as ds
import mosfet_deck as M
M.build()
D,S=M.DEV,M.SI
x,y=M.coords(S)
ds.node_solution(device=D,region=S,name='nidx'); ds.set_node_values(device=D,region=S,name='nidx',values=list(map(float,range(len(x)))))
ds.edge_from_node_model(device=D,region=S,node_model='nidx')
i0=np.array(ds.get_edge_model_values(device=D,region=S,name='nidx@n0')).astype(int)
i1=np.array(ds.get_edge_model_values(device=D,region=S,name='nidx@n1')).astype(int)
ec=np.array(ds.get_edge_model_values(device=D,region=S,name='EdgeCouple'))
L=2e-4
cuts=[0.1*L,0.25*L,0.5*L,0.75*L,0.9*L,0.98*L]
def cutcur():
    J=np.array(ds.get_edge_model_values(device=D,region=S,name='ElectronCurrent'))*ec
    out=[]
    for xc in cuts:
        inS=x<xc
        a,b=inS[i0],inS[i1]
        # flux leaving S (from n0 to n1 positive)
        I=np.sum(J[a&~b])-np.sum(J[~a&b])
        out.append(I)
    bodyI=ds.get_contact_current(device=D,contact='body',equation='ElectronContinuityEquation')
    return np.array(out), ds.get_circuit_node_value(solution='dcop',node='V1.I')
r={}
for v in (0.0,2.5e-4,-2.5e-4,1e-3,-1e-3):
    ds.circuit_alter(name='V1',value=v); M.solve_dc(); r[v]=cutcur(); print(v,r[v],flush=True)
for d in (2.5e-4,1e-3):
    print(d,'cut G',(r[d][0]-r[-d][0])/(2*d),'V1 G',(r[d][1]-r[-d][1])/(2*d))
