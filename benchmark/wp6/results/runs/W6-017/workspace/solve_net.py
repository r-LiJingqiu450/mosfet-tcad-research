import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
d=np.load('network.npz'); i0,i1,g=d['i0'],d['i1'],d['g']; N=len(d['psi'])
keep=g>0; i0,i1,g=i0[keep],i1[keep],g[keep]
print('g range',g.min(),g.max())
Lap=sp.coo_matrix((np.r_[-g,-g,g,g],(np.r_[i0,i1,i0,i1],np.r_[i1,i0,i0,i1])),shape=(N,N)).tocsr()
fixed=np.zeros(N,bool); val=np.zeros(N)
for k in ('src','bod'): fixed[d[k]]=True
fixed[d['drn']]=True; val[d['drn']]=1.0
assert not np.any(np.isin(d['drn'],np.r_[d['src'],d['bod']]))
fr=~fixed
A=Lap[fr][:,fr].tocsc(); b=-Lap[fr][:,fixed]@val[fixed]
# diagonal scaling
s=1/np.sqrt(A.diagonal()); As=sp.diags(s)@A@sp.diags(s)
lu=sla.splu(As.tocsc())
x=s*lu.solve(s*b)
# iterative refinement in longdouble
Al=A.astype(np.longdouble)
for it in range(5):
    r=(b.astype(np.longdouble)-Al@x.astype(np.longdouble))
    print('refine',it,float(np.max(np.abs(r))))
    x=(x.astype(np.longdouble)+s*lu.solve((s*r).astype(float))).astype(np.longdouble)
phi=val.astype(np.longdouble); phi[fr]=x
dphi=phi[i0]-phi[i1]
I_edge=g*dphi  # from n0 to n1
def flux(nodes):
    m=np.zeros(N,bool); m[nodes]=True
    out=np.sum(I_edge[m[i0]&~m[i1]])-np.sum(I_edge[m[i1]&~m[i0]])
    return float(out)
print('drain out-flux',flux(d['drn']))
print('source in',-flux(d['src']),'body in',-flux(d['bod']))
print('energy G',float(np.sum(g*dphi**2)))
# KCL residual at interior
