import pickle, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
d=pickle.load(open('jac.pkl','rb'))
JA,JB,FA,FB,xA,xB,E,C=[d[k] for k in 'JA JB FA FB xA xB E C'.split()]
iI=C['V1.I']; idb=C['drain_bias']
J=JA.tocsr()
r=1/np.array(abs(J).max(1).todense()).ravel()
J1=sp.diags(r)@J
c=1/np.array(abs(J1).max(0).todense()).ravel()
J2=(J1@sp.diags(c)).tocsc()
for opt in ['COLAMD','MMD_AT_PLUS_A']:
    try:
        lu=spl.splu(J2,permc_spec=opt); break
    except Exception as e: print(opt,e)
b=np.zeros(len(FA)); b[iI]=1.0
y=lu.solve(r*b); x=c*y
print('residual', np.max(np.abs(JA@x-b)))
print('G lin', -x[iI]/x[idb], x[idb])
print('G fd', -(xB[iI]-xA[iI])/(xB[idb]-xA[idb]))
np.save('xlin.npy',x)
