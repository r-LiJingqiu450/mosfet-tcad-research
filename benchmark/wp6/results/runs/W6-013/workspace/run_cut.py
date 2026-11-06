import sys, numpy as np
import mosfet_deck as M, devsim as ds
M.build()
D,R='mos','si'
for nm in ('x','y','node_index'):
    ds.edge_from_node_model(device=D, region=R, node_model=nm)
E=lambda n: np.array(ds.get_edge_model_values(device=D, region=R, name=n))
N=lambda n: np.array(ds.get_node_model_values(device=D, region=R, name=n))
x0,x1,i0,i1=E('x@n0'),E('x@n1'),E('node_index@n0').astype(int),E('node_index@n1').astype(int)
ec=E('EdgeCouple'); xn,yn=N('x'),N('y'); nn=len(xn)
xs=np.unique(xn); 
body=np.where(np.abs(yn-yn.max())<1e-15)[0]
drain=np.where((np.abs(xn-xn.max())<1e-15)&(yn<=M.MESH['contact_extent_cm']*(1+1e-12)))[0]
def currents():
    F=E('ElectronCurrent')*ec   # conventional current n0->n1, A/cm
    out=np.zeros(nn); np.add.at(out,i0,F); np.add.at(out,i1,-F)   # net conventional outflow from node into device
    res={}
    for xc in cuts:
        # current crossing x=xc from left to right
        m1=(x0<xc)&(x1>xc); m2=(x0>xc)&(x1<xc)
        Icut=F[m1].sum()-F[m2].sum()
        bR=body[xn[body]>xc]
        IbR=out[bR].sum()   # outflow from right body contact nodes into device
        # current into drain terminal from device side: I_D = -(outflow of drain nodes) ; via cut: drain+bodyR outflow = -Icut... 
        res[xc]=(Icut, IbR, -Icut - IbR)
    return res, out[drain].sum(), out[body].sum()
cuts=[0.5*(xs[k]+xs[k+1]) for k in np.searchsorted(xs,[0.2e-4,0.5e-4,0.8e-4])]
print('cuts',cuts)
vals={}
for v in [0.0, 1e-3,-1e-3,1e-2,-1e-2]:
    ds.circuit_alter(name='V1', value=v); M.solve_dc()
    r,od,ob=currents(); vals[v]=r
    print(v, {k:tuple('%.6e'%a for a in t) for k,t in r.items()}, 'drain_out %.4e body_out %.4e'%(od,ob),
          'V1.I %.6e'%ds.get_circuit_node_value(solution='dcop', node='V1.I'), flush=True)
for dv in (1e-3,1e-2):
    for xc in cuts:
        print('dv',dv,'xc',xc,'G_D(-(cut+bodyR) outflow)=%.6e'%(-(vals[dv][xc][2]-vals[-dv][xc][2])/(2*dv)) if False else
              'G=%.6e  Gcut=%.6e GbR=%.6e'%(((vals[dv][xc][2]-vals[-dv][xc][2])/(2*dv)),(vals[dv][xc][0]-vals[-dv][xc][0])/(2*dv),(vals[dv][xc][1]-vals[-dv][xc][1])/(2*dv)))
