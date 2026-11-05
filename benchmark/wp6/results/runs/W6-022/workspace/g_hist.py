import devsim as ds, mosfet_deck as M
M.build()
for dv in (1e-4, 2e-4):
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc()
    for h in r['history']: print(dv, h)
    print('I',ds.get_circuit_node_value(solution='dcop', node='V1.I'))
