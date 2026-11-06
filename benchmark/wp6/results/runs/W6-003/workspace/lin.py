import devsim as ds, mosfet_deck as M, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
M.build()
D,R='mos','si'
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
g=lambda n: np.array(ds.get_edge_model_values(device=D,region=R,name=n))
i0=g('node_index@n0').astype(int); i1=g('node_index@n1').astype(int)
EC=g('EdgeCouple'); EIL=g('EdgeInverseLength'); v=g('vdiff'); B=g('Bern01'); n0=g('Electrons@n0'); n1=g('Electrons@n1')
J=g('ElectronCurrent')
psi=np.array(ds.get_node_model_values(device=D,region=R,name='Potential'))
n=np.array(ds.get_node_model_values(device=D,region=R,name='Electrons'))
phin=psi-M.VT*np.log(n/M.NI)
print('max|phin| at equilibrium',np.abs(phin).max())
print('max |n1 e^v / n0 -1|', np.abs(n1*np.exp(v)/n0-1).max())
cond = M.Q*M.MU_N*EIL*EC*B*n0          # dJ/d(phin0-phin1) * EdgeCouple
print('edges',len(cond),'cond range',cond[cond>0].min(),cond.max(), 'zero-couple edges',(EC==0).sum())
N=len(psi)
def cnodes(c):
    s=set()
    for e in ds.get_element_node_list(device=D,region=R,contact=c):
        s.update(e)
    return np.array(sorted(s))
cd,cs,cb=cnodes('drain'),cnodes('source'),cnodes('body')
print(len(cd),len(cs),len(cb), set(cd)&set(cs), set(cd)&set(cb), set(cs)&set(cb))
L=sp.coo_matrix((np.concatenate([cond,cond,-cond,-cond]),(np.concatenate([i0,i1,i0,i1]),np.concatenate([i0,i1,i1,i0]))),shape=(N,N)).tocsr()
phi=np.zeros(N); phi[cd]=1.0
fixed=np.zeros(N,bool); fixed[cd]=fixed[cs]=fixed[cb]=True
fr=np.where(~fixed)[0]
Lff=L[fr][:,fr].tocsc(); rhs=-(L[fr][:,cd]@np.ones(len(cd)))
# diagonal scaling for conditioning
d=1/np.sqrt(Lff.diagonal()); S=sp.diags(d)
A=(S@Lff@S).tocsc(); lu=sla.splu(A)
y=lu.solve(d*rhs)
for k in range(5):
    r=d*rhs-A@y; y+=lu.solve(r); print('refine resid',np.abs(r).max())
phi[fr]=d*y
dphi=phi[i0]-phi[i1]
G_energy=np.sum(cond*dphi**2)
Id=(L@phi)[cd].sum(); Ib=(L@phi)[cb].sum(); Is=(L@phi)[cs].sum()
print('G_energy',G_energy,'I_drain',Id,'I_source',Is,'I_body',Ib)
np.savez('lin_result.npz',phi=phi,cond=cond,i0=i0,i1=i1)
