import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, pickle
import devsim as ds
import mosfet_deck as M
M.build()
D='mos'
vars_=[('si','Potential'),('si','Electrons'),('ox','Potential')]
cn=ds.get_circuit_node_list(); print('circuit nodes',cn, ds.get_circuit_solution_list())
def getx():
    return {k:np.array(ds.get_node_model_values(device=D,region=k[0],name=k[1])) for k in vars_}, {c:ds.get_circuit_node_value(solution='dcop',node=c) for c in cn}
def setx(s):
    for k,v in s[0].items(): ds.set_node_values(device=D,region=k[0],name=k[1],values=list(v))
    for c,v in s[1].items(): ds.set_circuit_node_value(solution='dcop',node=c,value=v)
def eqn():
    E={k:np.array(ds.get_equation_numbers(device=D,region=k[0],variable=k[1])) for k in vars_}
    C={c:ds.get_circuit_equation_number(node=c) for c in cn}
    return E,C
def vec(s,E,C,N):
    x=np.zeros(N)
    for k,v in s[0].items(): x[E[k]]=v
    for c,v in s[1].items(): x[C[c]]=v
    return x
A=getx()
dv=1e-5
ds.circuit_alter(name='V1',value=dv); M.solve_dc()
B=getx()
r=ds.get_matrix_and_rhs(format='csc'); print(r.keys())
E,C=eqn(); print(C)
def mat():
    r=ds.get_matrix_and_rhs(format='csc')
    st=r['static']
    if not hasattr(mat,'shown'):
        mat.shown=1; print(type(st), st.keys() if hasattr(st,'keys') else len(st)); print('perm', type(r['permutation']), len(r['permutation']), list(r['permutation'])[:10])
    n=len(st['ap'])-1
    J=sp.csc_matrix((np.array(st['av']),np.array(st['ai']),np.array(st['ap'])),shape=(n,n))
    return J, np.array(st['rhs']), r['permutation']
JB,FB,P=mat()
setx(A)
JA,FA,P2=mat()
print('perm', {k:(len(v), list(v)[:5]) for k,v in P.items()})
N=JA.shape[0]; print('N',N)
xA=vec(A,E,C,N); xB=vec(B,E,C,N); d=xB-xA
pickle.dump(dict(JA=JA,JB=JB,FA=FA,FB=FB,xA=xA,xB=xB,E=E,C=C,A=A,B=B),open('jac.pkl','wb'))
dF=FB-FA
for nm,Jd in (('JA',JA@d),('JB',JB@d),('avg',0.5*(JA+JB)@d)):
    print(nm,'sign check', np.dot(Jd,dF)/np.dot(Jd,Jd), 'max |dF-Jd|',np.max(np.abs(dF-Jd)), 'max|dF|',np.max(np.abs(dF)))
# linear response from JA
dx=spl.spsolve(JA.tocsc(), FA - (FA*0)) if False else None
