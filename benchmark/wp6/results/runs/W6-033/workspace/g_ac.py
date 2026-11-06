import mosfet_deck as M
import devsim as ds
M.build()
for f in (0.0, 1.0, 1e3):
    ds.solve(type='ac', frequency=f, solver_type='direct')
    print('f',f,'Y real', ds.get_circuit_node_value(solution='ssac_real', node='V1.I'),
          'imag', ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'), flush=True)
