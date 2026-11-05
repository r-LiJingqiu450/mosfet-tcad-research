import pickle, numpy as np
d=pickle.load(open('jac.pkl','rb'))
J=d['JA'].tocsr(); E=d['E']; C=d['C']
lab=np.empty(J.shape[0],object)
for k,v in E.items(): lab[v]=[f'{k[0]}:{k[1]}:{i}' for i in range(len(v))]
for c,v in C.items(): lab[v]=c
rn=np.array(abs(J).sum(1)).ravel(); cnn=np.array(abs(J).sum(0)).ravel()
print('zero rows', [lab[i] for i in np.where(rn==0)[0]][:20], 'zero cols', [lab[i] for i in np.where(cnn==0)[0]][:20])
for k,v in E.items(): print(k, len(v), v.min(), v.max(), len(np.unique(v)))
print(J[C['V1.I']], J[C['drain_bias']].nnz)
