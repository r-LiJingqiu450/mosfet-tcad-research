import numpy as np, devsim as ds, mosfet_deck as M, scipy.sparse as spm, scipy.sparse.linalg as sla
D,S=M.DEV,M.SI
M.build()
st=ds.get_matrix_and_rhs(format='csr')['static']; N=len(st['rhs'])
J=spm.csr_matrix((np.array(st['av']),np.array(st['ai']),np.array(st['ap'])),shape=(N,N))
eqp=np.array(ds.get_equation_numbers(device=D,region=S,variable='Potential'))
eqn=np.array(ds.get_equation_numbers(device=D,region=S,variable='Electrons'))
ii=ds.get_circuit_equation_number(node='V1.I')
b=np.zeros(N); b[ii]=1.0
lu=sla.splu(J.tocsc())
x=lu.solve(b)
for k in range(5):
    r=b-J@x; print('refine',k,'dI/dV',x[ii],'resid',np.abs(r).max()); x=x+lu.solve(r)
# row/col scaling solve
dr=1/np.sqrt(abs(J).max(axis=1).toarray().ravel()); dc=1/np.sqrt(abs(J).max(axis=0).toarray().ravel())
Js=spm.diags(dr)@J@spm.diags(dc); lu2=sla.splu(Js.tocsc()); y=lu2.solve(dr*b); x2=dc*y
for k in range(5):
    r=b-J@x2; print('scaled refine',k,x2[ii],np.abs(r).max()); x2=x2+dc*lu2.solve(dr*r)
for nm in ('x',):
    ds.edge_from_node_model(device=D, region=S, node_model=nm)
ev=lambda n: np.array(ds.get_edge_model_values(device=D, region=S, name=n))
x0,x1=ev('x@n0'),ev('x@n1'); C=ev('EdgeCouple')
# edge node indices
en=np.array(ds.get_element_node_list(device=D, region=S)) if False else None
import itertools
def node_idx():
    # use node_index model
    return np.array(ds.get_edge_model_values(device=D, region=S, name='node_index@n0')).astype(int), np.array(ds.get_edge_model_values(device=D, region=S, name='node_index@n1')).astype(int)
ds.edge_from_node_model(device=D, region=S, node_model='node_index')
i0,i1=node_idx()
for xx,lab in ((x,'lu'),(x2,'scaled')):
    dp=xx[eqp]; dn=xx[eqn]
    dJ=ev('ElectronCurrent:Potential@n0')*dp[i0]+ev('ElectronCurrent:Potential@n1')*dp[i1]+ev('ElectronCurrent:Electrons@n0')*dn[i0]+ev('ElectronCurrent:Electrons@n1')*dn[i1]
    out=[]
    for xc in [0.25e-4,1e-4,1.75e-4]:
        a=(x0<xc)&(x1>xc); bb=(x1<xc)&(x0>xc); out.append(np.sum(dJ[a]*C[a])-np.sum(dJ[bb]*C[bb]))
    print(lab,'dI/dV circuit',xx[ii],'cut fluxes',out)
