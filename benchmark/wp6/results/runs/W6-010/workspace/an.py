import pickle, numpy as np
d=pickle.load(open('jac.pkl','rb'))
JA,JB,FA,FB,xA,xB,E,C=[d[k] for k in 'JA JB FA FB xA xB E C'.split()]
dx=xB-xA; dF=FB-FA
J=0.5*(JA+JB)
res=dF + J@dx
scale=abs(J)@abs(dx)+np.abs(dF)+1e-300
rel=np.abs(res)/scale
lab=np.empty(len(dx),object)
for k,v in E.items():
    lab[v]=[f'{k[0]}:{k[1]}:{i}' for i in range(len(v))]
for c,v in C.items(): lab[v]=c
o=np.argsort(-rel)[:25]
for i in o: print(lab[i], 'res',res[i],'scale',scale[i],'rel',rel[i], 'dF',dF[i])
print('drain_bias row etc', [(c,res[v],dF[v],(J@dx)[v]) for c,v in C.items()])
