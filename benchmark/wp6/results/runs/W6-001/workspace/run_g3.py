import numpy as np
import mosfet_deck as M
import devsim as ds
M.build()
D,R=M.DEV,M.SI
for nm in ('PotentialNodeCharge:Electrons','PotentialNodeCharge:Potential','NCharge:Electrons'):
    print(nm, ds.get_node_model_values(device=D,region=R,name=nm)[:3])
for nm in ('ElectronCurrent:Electrons@n0','ElectronCurrent:Electrons@n1','ElectronCurrent:Potential@n0','ElectronCurrent:Potential@n1'):
    print(nm, ds.get_edge_model_values(device=D,region=R,name=nm)[:3])
print(ds.get_edge_model_values(device=D,region=R,name='ElectronCurrent')[:3])
ds.circuit_alter(name='V1',value=1e-3)
h=M.solve_dc(); print(h['history'])
for f in (1e-3,1.0,1e3):
    ds.solve(type='ac', frequency=f)
    print('ac',f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'))
