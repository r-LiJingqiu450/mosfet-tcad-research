import mosfet_deck as M, devsim as ds, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
M.build()
D,R='mos','si'
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
ds.edge_model(device=D,region=R,name='i0',equation='node_index@n0')
ds.edge_model(device=D,region=R,name='i1',equation='node_index@n1')
g=lambda n: np.array(ds.get_edge_model_values(device=D,region=R,name=n))
i0=g('i0').astype(int); i1=g('i1').astype(int)
EC=g('EdgeCouple'); EL=g('EdgeLength'); v=g('vdiff'); B=g('Bern01')
psi=np.array(ds.get_node_model_values(device=D,region=R,name='Potential'))
n=np.array(ds.get_node_model_values(device=D,region=R,name='Electrons'))
x,y=M.coords(R)
# check equilibrium: phi_n
phi=psi-M.VT*np.log(n/M.NI); print('max |phi_n| eq', np.abs(phi).max())
for nb_kind in ['psi','n']:
    n0 = M.NI*np.exp(psi[i0]/M.VT) if nb_kind=='psi' else n[i0]
    w = M.Q*M.MU_N*EC/EL*n0*B   # conductance per edge (S/cm)
    N=len(psi)
    A=sp.coo_matrix((np.concatenate([w,w,-w,-w]),(np.concatenate([i0,i1,i0,i1]),np.concatenate([i0,i1,i1,i0]))),shape=(N,N)).tocsr()
    nodes={c:np.array(ds.get_element_node_list(device=D,region=R,contact=c)).ravel() for c in ['source','drain','body']}
    nodes={c:np.unique(v_) for c,v_ in nodes.items()}
    print({c:len(v_) for c,v_ in nodes.items()})
    fixed=np.zeros(N,bool); val=np.zeros(N)
    for c in nodes: fixed[nodes[c]]=True
    val[nodes['drain']]=1.0
    free=~fixed
    Aff=A[free][:,free]; rhs=-A[free][:,fixed]@val[fixed]
    sol=val.copy(); sol[free]=spl.spsolve(Aff.tocsc(),rhs)
    r=A@sol
    print(nb_kind,'G drain',r[nodes['drain']].sum(),'source',r[nodes['source']].sum(),'body',r[nodes['body']].sum(), 'free resid', np.abs(r[free]).max())
# shifted: drain=0, source=body=-1  -> G_DD = -I_drain
val=np.zeros(N); val[nodes['source']]=-1; val[nodes['body']]=-1
sol=val.copy(); sol[free]=spl.spsolve(Aff.tocsc(),-A[free][:,fixed]@val[fixed])
# iterative refinement
for k in range(3):
    r=A@sol; sol[free]-=spl.spsolve(Aff.tocsc(), r[free])
r=A@sol
print('shifted: G_DD=',-r[nodes['drain']].sum(),'source',r[nodes['source']].sum(),'body',r[nodes['body']].sum())
