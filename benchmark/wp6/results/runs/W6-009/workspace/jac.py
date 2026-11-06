import numpy as np, mosfet_deck as M, devsim as ds
M.build()
for v in [1e-5]:
    ds.circuit_alter(name='V1',value=v)
    for k in range(6):
        e=M._errors(M._solve_one())
        print(k, {kk:(vv if not isinstance(vv,dict) else vv['abs']) for kk,vv in e.items()}, ds.get_circuit_node_value(node='V1.I',solution='dcop'))
