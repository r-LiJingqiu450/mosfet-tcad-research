import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
import mosfet_deck as M, devsim as ds
M.build()
r=ds.get_matrix_and_rhs(format='csc'); s=r['static']
row=np.array(r['permutation']['row']); cp=np.array(r['permutation']['copy'])
N=len(row); idx=np.arange(N)
print('row!=i', np.sum(row!=idx), 'copy unique', np.unique(cp)[:10], np.sum(cp!=0))
w=np.nonzero(row!=idx)[0]; print(w[:20], row[w[:20]], cp[w[:20]])
F=np.array(s['rhs']); print('F at top rows', np.argsort(-np.abs(F))[:10], np.sort(-np.abs(F))[:10])
for c in ('source','drain','body'):
    for eq in ('ElectronContinuityEquation','PotentialEquation'):
        print(c,eq, ds.get_contact_current(device='mos',contact=c,equation=eq))
# map si node indices for big residual
x,y=M.coords('si'); k=np.argmax(np.abs(F))-44160; print('node',k,x[k],y[k])
np.save('F0.npy',F)
