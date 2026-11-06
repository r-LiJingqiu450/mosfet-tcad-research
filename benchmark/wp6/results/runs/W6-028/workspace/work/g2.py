import numpy as np, math, scipy.sparse as sp, scipy.sparse.linalg as sla
import mosfet_deck as M, devsim as ds
M.build()
D,R=M.DEV,M.SI
g=lambda n: np.array(ds.get_node_model_values(device=D,region=R,name=n))
ge=lambda n: np.array(ds.get_edge_model_values(device=D,region=R,name=n))
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
i0=ge('node_index@n0').astype(int); i1=ge('node_index@n1').astype(int)
psi=g('Potential'); n=g('Electrons'); x,y=M.coords()
VT=M.VT; NI=M.NI
print('max |log(n/(ni e^psi/VT))|', np.max(np.abs(np.log(n)-np.log(NI)-psi/VT)))
ec=ge('EdgeCouple'); el=ge('EdgeLength')
a=psi[i0]/VT; b=psi[i1]/VT; v=a-b
# n_eff = ni * v e^a /(e^v-1) = ni*e^b * v/(1-e^-v)... compute stably: ni*exp(a)*B(v)
def B(v):
    out=np.empty_like(v); s=np.abs(v)<1e-8
    out[s]=1-v[s]/2; out[~s]=v[~s]/np.expm1(v[~s]); return out
neff=NI*np.exp(a)*B(v)
C=M.Q*M.MU_N*neff*ec/el
N=len(psi)
def contact_nodes(c):
    return np.array(ds.get_node_model_values(device=D,region=R,name='node_index'))[np.array(ds.get_node_model_values(device=D,region=R,name='node_index')).astype(int)] if False else None
# contact node sets
xs_min,xs_max=x.min(),x.max(); ce=M.MESH['contact_extent_cm']*(1+1e-12)
src=np.where((np.abs(x-xs_min)<1e-15)&(y<=ce))[0]; drn=np.where((np.abs(x-xs_max)<1e-15)&(y<=ce))[0]
body=np.where(np.abs(y-y.max())<1e-15)[0]
print(len(src),len(drn),len(body))
L=sp.coo_matrix((np.concatenate([C,C,-C,-C]),(np.concatenate([i0,i1,i0,i1]),np.concatenate([i0,i1,i1,i0]))),shape=(N,N)).tocsr()
fixed=np.zeros(N,bool); fixed[src]=fixed[drn]=fixed[body]=True
phi=np.zeros(N); phi[drn]=1.0
fr=~fixed
A=L[fr][:,fr]; rhs=-L[fr][:,fixed]@phi[fixed]
phi[fr]=sla.spsolve(A.tocsc(),rhs)
Iin=(L@phi)[drn].sum()
P=np.sum(C*(phi[i0]-phi[i1])**2)
Isrc=(L@phi)[src].sum(); Ib=(L@phi)[body].sum()
print('G drain flux',Iin,'G power',P,'src',Isrc,'body',Ib)
# mid cut
xm=0.5e-4; cut=((x[i0]<xm)&(x[i1]>xm))|((x[i1]<xm)&(x[i0]>xm))
sgn=np.where(x[i0]<xm,1,-1)
print('mid-cut flux', np.sum((C*(phi[i0]-phi[i1])*sgn)[cut]))
np.savez('work/net.npz',phi=phi,C=C)
