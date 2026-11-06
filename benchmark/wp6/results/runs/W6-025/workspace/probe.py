import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
import mosfet_deck as M, devsim as ds
M.build()
r=ds.get_matrix_and_rhs(format='csc')
print(type(r), r.keys() if isinstance(r,dict) else None)
for k,v in r.items():
    print(k, type(v), (v.keys() if isinstance(v,dict) else (len(v) if hasattr(v,'__len__') else v)))
for reg in ('si','ox'):
    for eq in ('PotentialEquation','ElectronContinuityEquation'):
        try:
            en=ds.get_equation_numbers(device='mos',region=reg,equation=eq); print(reg,eq,min(en),max(en),len(en))
        except Exception as e: print(reg,eq,e)
print(ds.get_circuit_node_list(), [ds.get_circuit_equation_number(node=n) for n in ds.get_circuit_node_list()] if hasattr(ds,'get_circuit_equation_number') else '')
