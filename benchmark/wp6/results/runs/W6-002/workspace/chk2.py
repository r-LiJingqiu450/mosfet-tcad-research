import devsim as ds, mosfet_deck as M
M.build()
d=1e-6
ds.circuit_alter(name='V1', value=d)
for k in range(5):
    M._solve_one()
    print(k, 'G=', -ds.get_circuit_node_value(solution='dcop', node='V1.I')/d, ds.get_circuit_node_value(solution='dcop', node='drain_bias'), flush=True)
