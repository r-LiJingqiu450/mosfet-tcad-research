import numpy as np
import mosfet_deck as M, devsim as ds
info=M.build()
xs,ys=M.coords()
def edges():
    x0=np.array(ds.get_edge_model_values(device='mos',region='si',name='x@n0')) if False else None
ds.edge_from_node_model(device='mos',region='si',node_model='x')
ds.edge_from_node_model(device='mos',region='si',node_model='y')
def cuts():
    J=np.array(ds.get_edge_model_values(device='mos',region='si',name='ElectronCurrent'))
    C=np.array(ds.get_edge_model_values(device='mos',region='si',name='EdgeCouple'))
    x0=np.array(ds.get_edge_model_values(device='mos',region='si',name='x@n0'))
    x1=np.array(ds.get_edge_model_values(device='mos',region='si',name='x@n1'))
    out=[]
    for xc in (-0.25e-4,0.0,0.5e-4,1.0e-4,1.25e-4):
        m=((x0<xc)&(x1>xc))|((x1<xc)&(x0>xc))
        s=np.sign(x1-x0)
        out.append(np.sum((J*C*s)[m]))
    return np.array(out)
print('cut0',cuts())
r={}
for dv in (1e-4,-1e-4):
    ds.circuit_alter(name='V1',value=dv); h=M.solve_dc(); print([f"{q['dpsi']:.1e}" for q in h['history']])
    r[dv]=(cuts(), ds.get_circuit_node_value(node='V1.I',solution='dcop'))
    print(dv, r[dv])
print('G cuts', (r[1e-4][0]-r[-1e-4][0])/2e-4, 'G V1', (r[1e-4][1]-r[-1e-4][1])/2e-4)
ds.circuit_alter(name='V1',value=0.); M.solve_dc()
for f in (1e-9,1e-6,1e-3,1.0):
    ds.solve(type='ac', frequency=f)
    print('ac',f, ds.get_circuit_node_value(node='V1.I',solution='ssac_real'), ds.get_circuit_node_value(node='V1.I',solution='ssac_imag'),flush=True)
