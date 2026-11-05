import mosfet_deck as M, devsim as ds, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
M.build()
D,R='mos','si'
g=lambda n: np.array(ds.get_edge_model_values(device=D,region=R,name=n))
nv=lambda n: np.array(ds.get_node_model_values(device=D,region=R,name=n))
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
i0=g('node_index@n0').astype(int); i1=g('node_index@n1').astype(int)
EC=g('EdgeCouple'); IL=g('EdgeInverseLength')
psi=nv('Potential'); n=nv('Electrons'); x=nv('x'); y=nv('y')
VT=M.VT; ni=M.NI
phin=psi-VT*np.log(n/ni)
print('phi_n range', phin.min(), phin.max(), 'n range', n.min(), n.max())
v=(psi[i0]-psi[i1])/VT
# log of n1*B(v)
with np.errstate(all='ignore'):
    logB=np.where(np.abs(v)<1e-8, -v/2, np.log(np.abs(v))-np.log(np.abs(np.expm1(v))))
lnn=np.log(ni)+psi/VT   # equilibrium n
geq=M.Q*M.MU_N*EC*IL*np.exp(lnn[i1]+logB)
gsol=M.Q*M.MU_N*EC*IL*np.exp(np.log(n[i1])+logB)
print('g range', geq.min(), geq.max(), 'max rel diff eq vs sol', np.max(np.abs(gsol/geq-1)[EC>0]))

N=len(psi)
ext=M.MESH['contact_extent_cm']*(1+1e-12)
src=(np.abs(x-x.min())<1e-15)&(y<=ext); drn=(np.abs(x-x.max())<1e-15)&(y<=ext); bod=np.abs(y-y.max())<1e-15
fixed=src|drn|bod
print('contact nodes', src.sum(), drn.sum(), bod.sum())
for label,G in (('eq',geq),('sol',gsol)):
    A=sp.coo_matrix((np.concatenate([G,G,-G,-G]),(np.concatenate([i0,i1,i0,i1]),np.concatenate([i0,i1,i1,i0]))),shape=(N,N)).tocsr()
    phi=np.zeros(N); phi[drn]=1.0
    fr=~fixed
    Aff=A[fr][:,fr]; b=-A[fr][:,fixed]@phi[fixed]
    phi[fr]=sla.spsolve(Aff.tocsc(),b)
    # residual refine
    for k in range(3):
        r=b-Aff@phi[fr]; phi[fr]+=sla.spsolve(Aff.tocsc(),r)
    J=G*(phi[i0]-phi[i1])  # current (in units of q-flux) from n0 to n1 per V
    # cut currents at several x
    out=[]
    for xc in (0.1e-4,0.3e-4,0.5e-4,0.7e-4,0.9e-4):
        a=x[i0]<xc; c=x[i1]<xc
        cr=a!=c
        out.append(np.sum(np.where(a[cr],J[cr],-J[cr])))
    # contact currents: sum of flux out of fixed nodes into... 
    Id=sum(np.sum(J[drn[i0]&~drn[i1]])-np.sum(J[drn[i1]&~drn[i0]]) for _ in [0])
    Ib=np.sum(J[bod[i1]&~bod[i0]])-np.sum(J[bod[i0]&~bod[i1]])
    Is=np.sum(J[src[i1]&~src[i0]])-np.sum(J[src[i0]&~src[i1]])
    print(label,'cut fluxes (source<-drain, S/cm)',out,'drain out',Id,'body in',Ib,'source in',Is)
k=np.argmin(np.abs(x-0.5e-4)+np.abs(y))
print('mid surface psi', psi[k], 'band bending', psi[k]-M.PSI_B0, 'n_s', n[k])
col=np.abs(x-x[k])<1e-15
o=np.argsort(y[col]); print('N_inv approx', np.trapezoid(n[col][o],y[col][o]))
G=geq
A=sp.coo_matrix((np.concatenate([G,G,-G,-G]),(np.concatenate([i0,i1,i0,i1]),np.concatenate([i0,i1,i1,i0]))),shape=(N,N)).tocsr()
phi=np.zeros(N); phi[drn]=1.0; fr=~fixed
phi[fr]=sla.spsolve(A[fr][:,fr].tocsc(),-A[fr][:,fixed]@phi[fixed])
J=G*(phi[i0]-phi[i1])
xc=0.5e-4; a=x[i0]<xc; c=x[i1]<xc; cr=a!=c
Jc=np.where(a[cr],J[cr],-J[cr]); yc=0.5*(y[i0]+y[i1])[cr]
o=np.argsort(yc); cum=np.cumsum(Jc[o])
for f in (0.01,0.1,0.5,0.9,0.99): print(f, yc[o][np.searchsorted(cum/cum[-1],f)])
print('edges crossing', cr.sum(), 'top contributions', list(zip(yc[o][:5],Jc[o][:5])))
# check diagonal couples
diag=(x[i0]!=x[i1])&(y[i0]!=y[i1]); print('diag couple max', np.abs(EC[diag]).max(), 'neg couples', (EC<0).sum())
