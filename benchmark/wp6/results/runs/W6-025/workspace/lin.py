import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
import mosfet_deck as M, devsim as ds
M.build()
def get():
    r=ds.get_matrix_and_rhs(format='csc'); s=r['static']
    N=len(s['rhs'])
    J=sp.csc_matrix((s['av'],s['ai'],s['ap']),shape=(N,N))
    return J, np.array(s['rhs']), r
J0,F0,r=get()
print('perm', {k:(len(v), v[:5]) for k,v in r['permutation'].items()})
d=1e-4
ds.circuit_alter(name='V1', value=d)
J1,F1,_=get()
dF=(F1-F0)/d
nz=np.nonzero(dF)[0]; print('dF nonzero rows', nz, dF[nz])
print('|F0| max', np.abs(F0).max(), 'argmax', np.abs(F0).argmax(), '|F1-F0| max', np.abs(F1-F0).max())
lu=sla.splu(J0.tocsc())
dx=-lu.solve(dF)
print('dI/dV (linear) V1.I row', dx[85376], 'drain_bias', dx[85377])
# Newton step check: compare with DC solve update
