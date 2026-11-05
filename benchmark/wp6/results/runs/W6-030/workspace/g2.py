import time, numpy as np, devsim as ds, mosfet_deck as M
info=M.build(); 
print('I0', ds.get_circuit_node_value(solution='dcop', node='V1.I'))
ds.circuit_alter(name='V1', value=1e-2); r=M.solve_dc()
for h in r['history']: print(h)
print('I', ds.get_circuit_node_value(solution='dcop', node='V1.I'))
