import mosfet_deck as M, devsim as ds
M.build()
print('dcop V1.I', ds.get_circuit_node_value(solution='dcop', node='V1.I'))
for f in [1e-3, 1.0, 1e3]:
    ds.solve(type='ac', frequency=f)
    print('f',f,'Re I',ds.get_circuit_node_value(solution='ssac_real', node='V1.I'),'Im I',ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
