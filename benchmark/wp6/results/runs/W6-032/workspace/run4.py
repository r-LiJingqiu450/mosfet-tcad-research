import numpy as np
import mosfet_deck as M, devsim as ds
M.build()
for k in range(3):
    M._solve_one(); ds.solve(type='ac', frequency=1e-6)
    print('AC', ds.get_circuit_node_value(node='V1.I',solution='ssac_real'),flush=True)
r={}
for dv in (1e-5,-1e-5,3e-4,-3e-4):
    ds.circuit_alter(name='V1',value=dv); M.solve_dc(); r[dv]=ds.get_circuit_node_value(node='V1.I',solution='dcop')
for d in (1e-5,3e-4): print('FD',d,-(r[d]-r[-d])/(2*d))
