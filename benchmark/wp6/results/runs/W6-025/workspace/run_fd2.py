import mosfet_deck as M, devsim as ds
M.build()
res={}
for v in [1e-5,-1e-5,3e-4,-3e-4]:
    ds.circuit_alter(name='V1', value=v); M.solve_dc()
    res[v]=ds.get_circuit_node_value(node='V1.I', solution='dcop')
for a in (1e-5,3e-4):
    print('step',a,'G',-(res[a]-res[-a])/(2*a))
