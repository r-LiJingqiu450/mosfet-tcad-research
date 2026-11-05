import devsim as ds, mosfet_deck as M
M.build()
for f in (0.0, 1.0):
    ds.solve(type='ac', frequency=f)
    print('f',f,'Re',ds.get_circuit_node_value(solution='ssac_real', node='V1.I'),'Im',ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
