import pickle, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
d0=pickle.load(open('work/mat0.pkl','rb')); d1=pickle.load(open('work/mat1.pkl','rb'))
s=d0['static']; N=len(s['rhs'])
J=sp.csc_matrix((s['av'],s['ai'],s['ap']),shape=(N,N))
r0=np.array(s['rhs']); r1=np.array(d1['static']['rhs'])
print('perm row', d0['permutation']['row'][:5], 'nontrivial', np.sum(np.array(d0['permutation']['row'])!=np.arange(N)), np.sum(np.array(d0['permutation']['copy'])!=-1) if True else '')
db=r1-r0; nz=np.flatnonzero(db); print('db nz', nz, db[nz]); print('|r0|max', np.abs(r0).max(), 'at', np.argmax(np.abs(r0)))
lu=spl.splu(J, permc_spec='COLAMD')
for rhs,name in ((db,'db'),(r1,'r1')):
    x=lu.solve(rhs)
    print(name,'dx V1.I', x[273312], 'drain_bias', x[273313])
    # iterative refinement with longdouble residual
    Jl=sp.csr_matrix(J)
    for it in range(4):
        res = rhs.astype(np.longdouble) - (Jl.astype(np.longdouble) @ x.astype(np.longdouble))
        print('  it',it,'|res|',float(np.abs(res).max()), 'V1.I', x[273312])
        x = x + lu.solve(np.asarray(res, dtype=float))
