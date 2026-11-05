import numpy as np, devsim as ds, mosfet_deck as M, scipy.sparse.linalg as sla
D=M.DEV
M.build()
def state():
    return dict(psi_si=np.array(ds.get_node_model_values(device=D,region='si',name='Potential')),
                n_si=np.array(ds.get_node_model_values(device=D,region='si',name='Electrons')),
                psi_ox=np.array(ds.get_node_model_values(device=D,region='ox',name='Potential')),
                vd=ds.get_circuit_node_value(node='drain_bias',solution='dcop'),
                i=ds.get_circuit_node_value(node='V1.I',solution='dcop'))
s0=state()
mats=ds.get_matrix_and_rhs(format='csr'); print(mats.keys())
print({k:type(v) for k,v in mats.items()}, mats['static'].keys() if isinstance(mats['static'],dict) else None); import scipy.sparse as spm; st=mats['static']; n=len(st['rhs']); J=spm.csr_matrix((np.array(st['av']),np.array(st['ai']),np.array(st['ap'])),shape=(n,n)); print('perm', {k:(type(v), len(v) if hasattr(v,'__len__') else v) for k,v in mats['permutation'].items()}); r0=np.array(mats['static']['rhs']); print(J.shape, np.abs(r0).max())
N=J.shape[0]
eq={}
eq['psi_si']=np.array(ds.get_equation_numbers(device=D,region='si',variable='Potential'))
eq['n_si']=np.array(ds.get_equation_numbers(device=D,region='si',variable='Electrons'))
eq['psi_ox']=np.array(ds.get_equation_numbers(device=D,region='ox',variable='Potential'))
eq['vd']=ds.get_circuit_equation_number(node='drain_bias'); eq['i']=ds.get_circuit_equation_number(node='V1.I')
print({k:(np.min(v),np.max(v)) for k,v in eq.items()})
d=1e-4
ds.circuit_alter(name='V1',value=d)
r1=np.array(ds.get_matrix_and_rhs(format='csr')['static']['rhs']); b=(r1-r0)/d
print('b nonzero', np.nonzero(b)[0], b[np.nonzero(b)[0]])
M.solve_dc(); sp=state()
ds.circuit_alter(name='V1',value=-d); M.solve_dc(); sm=state()
def vec(sa,sb):
    x=np.zeros(N)
    for k in eq:
        x[eq[k]]=(np.asarray(sa[k])-np.asarray(sb[k]))/(2*d) if k not in('psi_ox',) else 0
    # oxide: equation numbers may overlap interface; assign anyway
    x[eq['psi_ox']]=(sa['psi_ox']-sb['psi_ox'])/(2*d)
    return x
dx=vec(sp,sm)
print('FD dI/dV', dx[eq['i']], 'dVd', dx[eq['vd']])
res=J@dx+b
print('|J dx + b| max', np.abs(res).max(), ' |J dx - b| max', np.abs(J@dx-b).max())
xl=sla.spsolve(J.tocsc(), -b)
print('linear: dI/dV', xl[eq['i']], 'dVd', xl[eq['vd']])
row=np.array(mats['permutation']['row']); cp=np.array(mats['permutation']['copy'])
idx=np.nonzero(row!=np.arange(N))[0]; print('nperm', len(idx), idx[:5], row[idx[:5]], cp[idx[:5]]); print('copy uniq', np.unique(cp, return_counts=True))
print('row uniq count', len(np.unique(row)))
keep=row!=-1
for sgn in (1,-1):
    R=J@dx+sgn*b
    print('sgn',sgn,'max keep',np.abs(R[keep]).max(),'max drop',np.abs(R[~keep]).max())
    big=np.argsort(-np.abs(R*keep))[:8]; print(big, R[big], row[big])
print('r0 keep max', np.abs(r0[keep]).max(), 'drop', np.abs(r0[~keep]).max())
# apply permutation: rows with row=-1 removed; is row[i]==i for others?
print('nontrivial dest', np.sum((row!=-1)&(row!=np.arange(N))))
