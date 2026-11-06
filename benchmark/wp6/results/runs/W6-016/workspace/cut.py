import numpy as np
import mosfet_deck as M
import devsim as ds
D,R='mos','si'
M.build()
for v in ('x','y'): ds.edge_from_node_model(device=D,region=R,node_model=v)
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
def ev(n): return np.array(ds.get_edge_model_values(device=D,region=R,name=n))
x0,x1,y0,y1=ev('x@n0'),ev('x@n1'),ev('y@n0'),ev('y@n1')
i0,i1=ev('node_index@n0').astype(int),ev('node_index@n1').astype(int)
ec=ev('EdgeCouple')
xs=np.array(M.MESH['x_cm']); ymax=1.2e-4
def currents():
    J=ev('ElectronCurrent')*ec   # conventional current n0->n1
    out={}
    for xc in (0.2e-4,0.5e-4,0.8e-4):
        k=np.searchsorted(xs,xc); xm=0.5*(xs[k-1]+xs[k])
        a=(x0<xm)&(x1>xm); b=(x0>xm)&(x1<xm)
        out[f'cut{xc*1e4:.1f}']=J[a].sum()-J[b].sum()   # current flowing +x (drain->source is -x)
    # body contact current: edges with exactly one node at y=ymax
    on0=np.abs(y0-ymax)<1e-15; on1=np.abs(y1-ymax)<1e-15
    # current into body nodes from interior
    a=on1&~on0; b=on0&~on1
    xb=np.where(a,x1,x0)
    for xm in (0.5e-4,):
        out['body_right']=J[a&(xb>xm)].sum()-J[b&(xb>xm)].sum()
        out['body_left']=J[a&(xb<=xm)].sum()-J[b&(xb<=xm)].sum()
    # drain contact nodes
    xmax=xs[-1]; cy=M.MESH['contact_extent_cm']*(1+1e-12)
    d0=(np.abs(x0-xmax)<1e-15)&(y0<=cy); d1=(np.abs(x1-xmax)<1e-15)&(y1<=cy)
    out['drain_out']=J[d0&~d1].sum()-J[d1&~d0].sum()   # current leaving drain nodes into device
    out['V1.I']=ds.get_circuit_node_value(solution='dcop',node='V1.I')
    return out
res={}
for v in (1e-3,-1e-3,0.0):
    ds.circuit_alter(name='V1',value=v); M.solve_dc(); res[v]=currents(); print(v,res[v],flush=True)
for k in res[0.0]:
    print(k,'G=',(res[1e-3][k]-res[-1e-3][k])/2e-3)
