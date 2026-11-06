import sys, time; sys.path.insert(0,'.')
import numpy as np, scipy.sparse as sp
import devsim as ds, mosfet_deck as M
M.build()
d = ds.get_matrix_and_rhs(format='csc')
def show(x, pre=''):
    if isinstance(x, dict):
        for k,v in x.items(): print(pre+str(k), type(v).__name__, (len(v) if hasattr(v,'__len__') and not isinstance(v,str) else v)); show(v, pre+'  ') if isinstance(v,dict) else None
show(d)
import pickle; pickle.dump(d, open('work/mat0.pkl','wb'))
ds.circuit_alter(name='V1', value=1e-5)
d1 = ds.get_matrix_and_rhs(format='csc'); pickle.dump(d1, open('work/mat1.pkl','wb'))
for nm in ('V1.I','drain_bias'):
    try: print(nm, ds.get_circuit_equation_number(node=nm))
    except Exception as e: print(nm, e)
