import mosfet_deck as M, devsim as ds
M.build()
I=lambda: ds.get_circuit_node_value(solution='dcop', node='V1.I')
I0=I(); dv=1e-5
ds.circuit_alter(name='V1', value=dv)
for k in range(4):
    M._solve_one(); print('iter',k+1,'G',(I()-I0)/dv, flush=True)
