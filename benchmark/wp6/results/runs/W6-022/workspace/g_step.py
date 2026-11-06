import devsim as ds, mosfet_deck as M
M.build()
dv=1e-6
ds.circuit_alter(name='V1', value=dv)
for k in range(4):
    M._solve_one(); print(k, ds.get_circuit_node_value(solution='dcop', node='V1.I')/dv, flush=True)
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
ds.solve(type='ac', frequency=0.0)
print('ac', ds.get_circuit_node_value(solution='ssac_real', node='V1.I'))
ds.circuit_alter(name='V1', value=dv)
for k in range(3):
    M._solve_one(); print(k, ds.get_circuit_node_value(solution='dcop', node='V1.I')/dv, flush=True)
