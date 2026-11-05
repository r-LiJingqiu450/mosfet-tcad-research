import devsim as ds, mosfet_deck as M
M.build()
for v in (1e-3, 2e-3):
    ds.circuit_alter(name='V1', value=v); r=M.solve_dc()
    for h in r['history']: print(v, h)
    print('I', ds.get_circuit_node_value(node='V1.I', solution='dcop'))
