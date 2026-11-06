import devsim as ds, mosfet_deck as M
M.build()
for dv in [1e-6,-1e-6]:
    ds.circuit_alter(name='V1', value=dv); M._solve_one()
    print('ONE', dv, ds.get_circuit_node_value(solution='dcop', node='V1.I')/dv)
    M._solve_one(); print('TWO', dv, ds.get_circuit_node_value(solution='dcop', node='V1.I')/dv)
    ds.circuit_alter(name='V1', value=0.); M.solve_dc()
for f in [1e-6, 1e-3]:
    ds.solve(type='ac', frequency=f)
    print('AC', f, ds.get_circuit_node_value(solution='ssac_real', node='V1.I'))
