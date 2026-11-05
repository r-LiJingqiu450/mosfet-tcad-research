import numpy as np
import mosfet_deck as M
import devsim as ds
D,R='mos','si'
M.build()
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
ds.edge_from_node_model(device=D,region=R,node_model='x')
g=lambda n: np.array(ds.get_edge_model_values(device=D,region=R,name=n))
i0=g('node_index@n0').astype(int); i1=g('node_index@n1').astype(int)
EC=g('EdgeCouple'); EL=g('EdgeLength')
x=np.array(ds.get_node_model_values(device=D,region=R,name='x'))
y=np.array(ds.get_node_model_values(device=D,region=R,name='y'))
ybot=y.max(); body=np.abs(y-ybot)<1e-15
xs=np.unique(x)
cuts=[0.5*(xs[k]+xs[k+1]) for k in np.searchsorted(xs,[0.1e-4,0.3e-4,0.5e-4,0.7e-4,0.9e-4])]
LD=np.longdouble
def bern(v):
    v=v.astype(LD); out=np.empty_like(v)
    small=np.abs(v)<1e-6
    out[small]=1-v[small]/2+v[small]**2/12
    vs=v[~small]; out[~small]=vs/np.expm1(vs)
    return out
def currents():
    psi=np.array(ds.get_node_model_values(device=D,region=R,name='Potential'),dtype=LD)
    n=np.array(ds.get_node_model_values(device=D,region=R,name='Electrons'),dtype=LD)
    v=(psi[i0]-psi[i1])/LD(M.VT)
    Bv=bern(v); Bm=bern(-v)
    # current n0->n1 conventional per deck: q mu VT/L*(n1 B(-v) - n0 B(v)) * EdgeCouple
    J=LD(M.Q*M.MU_N*M.VT)/EL.astype(LD)*(n[i1]*Bm-n[i0]*Bv)*EC.astype(LD)
    # node 'outflow' sum like residual: edge contributes +J to n0, -J to n1
    res=np.zeros(len(x),dtype=LD); np.add.at(res,i0,J); np.add.at(res,i1,-J)
    out={}
    for xc in cuts:
        Dn=x>xc
        fl=np.sum(np.where(Dn[i0]&~Dn[i1],J,0))-np.sum(np.where(Dn[i1]&~Dn[i0],J,0))
        Ib=np.sum(res[body&Dn])
        out[xc]=(float(fl),float(Ib),float(fl-Ib))
    return out
V1=lambda: ds.get_circuit_node_value(node='V1.I',solution='dcop')
r0=currents(); print('0',V1(),r0)
res={}
for dv in [1e-4,-1e-4,1e-3,-1e-3]:
    ds.circuit_alter(name='V1',value=dv); M.solve_dc()
    res[dv]=(V1(),currents()); print(dv,res[dv],flush=True)
for h in [1e-4,1e-3]:
    print('h',h,'G_V1',-(res[h][0]-res[-h][0])/(2*h))
    for xc in cuts:
        a=res[h][1][xc]; b=res[-h][1][xc]
        print('  cut %.3e'%xc,[ (a[k]-b[k])/(2*h) for k in range(3)])
