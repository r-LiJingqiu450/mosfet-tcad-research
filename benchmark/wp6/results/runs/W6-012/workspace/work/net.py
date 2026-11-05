import sys, math; sys.path.insert(0,'.')
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import mosfet_deck as M, devsim as ds
M.build()
D,R='mos','si'
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
i0=np.array(ds.get_edge_model_values(device=D,region=R,name='node_index@n0')).astype(int)
i1=np.array(ds.get_edge_model_values(device=D,region=R,name='node_index@n1')).astype(int)
ec=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
il=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeInverseLength'))
psi=np.array(ds.get_node_model_values(device=D,region=R,name='Potential'))
n=np.array(ds.get_node_model_values(device=D,region=R,name='Electrons'))
VT=M.VT; NI=M.NI
neq=NI*np.exp(psi/VT)
print('max rel |n/neq-1|',np.max(np.abs(n/neq-1)))
v=(psi[i0]-psi[i1])/VT
def B(x):
    x=np.asarray(x,float); out=np.empty_like(x); s=np.abs(x)<1e-8
    out[s]=1-x[s]/2; out[~s]=x[~s]/np.expm1(x[~s]); return out
# g = q mu/h * couple * B(v) * n0, use log form for precision
for label,nn in (('neq',neq),('n',n)):
  g=M.Q*M.MU_N*il*ec*B(v)*nn[i0]
  g1=M.Q*M.MU_N*il*ec*B(-v)*nn[i1]
  print(label,'asym max', np.max(np.abs(g-g1)/np.maximum(np.abs(g),1e-300)))
  g=0.5*(g+g1)
  N=len(psi)
  L=sp.coo_matrix((np.concatenate([g,g,-g,-g]),(np.concatenate([i0,i1,i0,i1]),np.concatenate([i0,i1,i1,i0]))),shape=(N,N)).tocsr()
  def cnodes(c):
    return np.array(ds.get_node_model_values(device=D,region=R,name='node_index'))[[]] if False else None
  # contact nodes
  x,y=M.coords(R)
  ce=M.MESH['contact_extent_cm']*(1+1e-12)
  src=(np.abs(x-x.min())<1e-15)&(y<=ce); drn=(np.abs(x-x.max())<1e-15)&(y<=ce); body=np.abs(y-y.max())<1e-15
  print(src.sum(),drn.sum(),body.sum())
  fixed=src|drn|body; free=~fixed
  phi=np.zeros(N); phi[drn]=1.0
  A=L[free][:,free]; b=-L[free][:,fixed]@phi[fixed]
  # diagonal scaling
  d=np.sqrt(A.diagonal()); S=sp.diags(1/d)
  z=spl.spsolve((S@A@S).tocsc(), S@b); phi[free]=S@z
  res=A@phi[free]-b; print('rel res',np.linalg.norm(res)/np.linalg.norm(b))
  Inj=L@phi
  print(label,'G drain inj',Inj[drn].sum(),'src',Inj[src].sum(),'body',Inj[body].sum(),'energy',np.sum(g*(phi[i0]-phi[i1])**2))
np.savez('work/net.npz',i0=i0,i1=i1,ec=ec,il=il,psi=psi,n=n,x=x,y=y)
