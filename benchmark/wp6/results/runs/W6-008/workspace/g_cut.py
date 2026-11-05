import sys, numpy as np
import mosfet_deck as M, devsim as ds
M.build()
D,R=M.DEV,M.SI
ds.edge_from_node_model(device=D, region=R, node_model='x')
xs=np.array(M.MESH['x_cm'])
X0=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n0'))
X1=np.array(ds.get_edge_model_values(device=D,region=R,name='x@n1'))
EC=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
cuts=[xs[k] for k in (20,79,150,218,290,358,420)]
def currents():
    J=np.array(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent'))
    out=[]
    for xc in cuts:
        xc2=xc+1e-12
        m01=(X0<xc2)&(X1>xc2); m10=(X1<xc2)&(X0>xc2)
        out.append(np.sum(J[m01]*EC[m01])-np.sum(J[m10]*EC[m10]))
    return np.array(out), ds.get_circuit_node_value(solution='dcop', node='V1.I'), ds.get_contact_current(device=D,contact='drain',equation='ElectronContinuityEquation')
print('cuts um', np.array(cuts)*1e4)
print('V=0', currents())
for dv in [float(a) for a in sys.argv[1:]]:
    r={}
    for s in (1,-1):
        ds.circuit_alter(name='V1', value=s*dv); M.solve_dc(); r[s]=currents()
        print('V',s*dv, r[s], flush=True)
    print('dv',dv,'Gcut',(r[1][0]-r[-1][0])/(2*dv),'G_V1',(r[1][1]-r[-1][1])/(2*dv),'G_cont',(r[1][2]-r[-1][2])/(2*dv), flush=True)
