import devsim as ds, mosfet_deck as M
M.build()
ds.circuit_alter(name='V1', value=1e-3); r=M.solve_dc()
for h in r['history']: print(h)
print(ds.get_circuit_node_value(node='V1.I'))
