import mosfet_deck as M, devsim as ds
M.build(); I={}
for dv in [1e-5,-1e-5,2e-4,-2e-4]:
    ds.circuit_alter(name='V1',value=dv); M.solve_dc(); I[dv]=ds.get_circuit_node_value(node='V1.I',solution='dcop')
for h in [1e-5,2e-4]: print('h',h,'G',-(I[h]-I[-h])/(2*h))
