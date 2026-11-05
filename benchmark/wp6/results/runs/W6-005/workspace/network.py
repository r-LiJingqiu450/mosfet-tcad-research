import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, devsim as ds, mosfet_deck as M
M.build()
D,R='mos','si'
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
g=lambda n: np.array(ds.get_edge_model_values(device=D,region=R,name=n))
i0=g('node_index@n0').astype(int); i1=g('node_index@n1').astype(int)
EC=g('EdgeCouple'); IL=g('EdgeInverseLength'); B=g('Bern01'); v=g('vdiff')
n=np.array(ds.get_node_model_values(device=D,region=R,name='Electrons'))
psi=np.array(ds.get_node_model_values(device=D,region=R,name='Potential'))
phi=psi-M.VT*np.log(n/M.NI); print('phi_n range', phi.min(), phi.max())
# SG: J = q mu Vt IL (n1 B(-v) - n0 B(v)); dJ/dphi0 = q mu IL n0 B(v), linear in dphi
ge=M.Q*M.MU_N*IL*EC*n[i0]*B
ge2=M.Q*M.MU_N*IL*EC*n[i1]*np.exp(v)*B  # equal at equilibrium
print('ge consistency', np.max(np.abs(ge-ge2)/np.maximum(ge,1e-300)))
N=len(n)
Gm=sp.coo_matrix((np.r_[ge,ge,-ge,-ge],(np.r_[i0,i1,i0,i1],np.r_[i0,i1,i1,i0])),shape=(N,N)).tocsr()
cn=lambda c: np.array(sorted(set(ds.get_element_node_list(device=D,region=R,contact=c) and [k for e in ds.get_element_node_list(device=D,region=R,contact=c) for k in e])))
dr,so,bo=cn('drain'),cn('source'),cn('body')
print(len(dr),len(so),len(bo))
fixed=np.zeros(N,bool); fixed[dr]=fixed[so]=fixed[bo]=True
x=np.zeros(N); x[dr]=1.0
free=~fixed
A=Gm[free][:,free]; b=-Gm[free][:,fixed]@x[fixed]
x[free]=spl.spsolve(A.tocsc(),b)
Iin=Gm@x
print('G_DD network', Iin[dr].sum(), ' source', Iin[so].sum(), ' body', Iin[bo].sum())
xs,ys=M.coords(R)
print('residual max free', np.abs(Iin[free]).max(), 'ge range', ge.min(), ge.max())
E=np.sum(ge*(x[i0]-x[i1])**2); print('G energy', E)
xm=1e-4
for xm in (0.5e-4,1e-4,1.5e-4):
    cross=(xs[i0]<xm)!=(xs[i1]<xm)
    s=np.where(xs[i0]>=xm,1,-1)  # current from right to left
    Icut=np.sum((ge*(x[i0]-x[i1])*s)[cross])  # flow from n0 to n1 = ge*(x0-x1)
    # body current on right side flows drain->body
    bR=bo[xs[bo]>=xm]
    print('cut',xm,Icut,'bodyR',-Iin[bR].sum(),'total',Icut-Iin[bR].sum())
# iterative refinement
for k in range(3):
    r=np.zeros(N); r=-(Gm@x); r[fixed]=0
    print('refine',k,np.abs(r).max())
    x[free]+=spl.spsolve(A.tocsc(),r[free])
Iin=Gm@x
print('after refine: drain',Iin[dr].sum(),'source',Iin[so].sum(),'body',Iin[bo].sum(),'energy',np.sum(ge*(x[i0]-x[i1])**2))
np.savez('network_sol.npz',x=x,xs=xs,ys=ys,ge=ge,i0=i0,i1=i1,n=n,psi=psi)
