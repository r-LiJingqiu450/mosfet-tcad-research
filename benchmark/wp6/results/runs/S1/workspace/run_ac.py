import json, numpy as np
import mosfet_deck as M, devsim as ds
m=json.load(open('mesh.json')); xs=np.array(m['x_cm']); ys=np.array(m['y_si_cm'])
print('gate cols x:', xs[56], xs[204], 'y_si[43]', ys[43], m['contact_extent_cm'])
M.build()
for f in [1e-6, 1.0]:
    ds.solve(type='ac', frequency=f)
    print('AC f',f, ds.get_circuit_node_value(node='V1.I',solution='ssac_real'), ds.get_circuit_node_value(node='V1.I',solution='ssac_imag'))
