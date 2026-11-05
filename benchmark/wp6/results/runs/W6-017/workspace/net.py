import devsim as ds, mosfet_deck as M, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
M.build()
D,R='mos','si'
psi=np.array(ds.get_node_model_values(device=D,region=R,name='Potential'))
nsol=np.array(ds.get_node_model_values(device=D,region=R,name='Electrons'))
VT=M.VT
print('max |log(n/(ni exp(psi/VT)))|', np.max(np.abs(np.log(nsol/M.NI)-psi/VT)))
en=np.array(ds.get_element_node_list(device=D,region=R)) if False else None
el=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeLength'))
ec=np.array(ds.get_edge_model_values(device=D,region=R,name='EdgeCouple'))
ds.edge_from_node_model(device=D,region=R,node_model='node_index')
i0=np.rint(ds.get_edge_model_values(device=D,region=R,name='node_index@n0')).astype(int)
i1=np.rint(ds.get_edge_model_values(device=D,region=R,name='node_index@n1')).astype(int)
v=(psi[i0]-psi[i1])/VT
def logB(v):
    out=np.empty_like(v); s=np.abs(v)<1e-6
    out[s]=-v[s]/2+v[s]**2/24
    p=(~s)&(v>0); out[p]=np.log(v[p])-v[p]-np.log1p(-np.exp(-v[p]))
    m=(~s)&(v<0); out[m]=np.log(-v[m])-np.log1p(-np.exp(v[m]))
    return out
logg=np.log(M.Q*M.MU_N*M.NI)+np.log(np.maximum(ec,1e-300))+(-np.log(el))+psi[i0]/VT+logB(v)
g=np.where(ec>0,np.exp(logg),0.0)
# symmetry check
logg2=np.log(M.Q*M.MU_N*M.NI)+psi[i1]/VT+logB(-v)
print('sym check',np.max(np.abs((logg2-(logg-np.log(np.maximum(ec,1e-300))+np.log(el)))[ec>0])))
N=len(psi)
def cnodes(c): return np.array(ds.get_element_node_list(device=D,region=R,contact=c)).ravel()
src=np.unique(cnodes('source')); drn=np.unique(cnodes('drain')); bod=np.unique(cnodes('body'))
print(len(src),len(drn),len(bod), 'edges',len(g),'zero couple',np.sum(ec<=0))
np.savez('network.npz',i0=i0,i1=i1,g=g,src=src,drn=drn,bod=bod,psi=psi,ec=ec,el=el)
