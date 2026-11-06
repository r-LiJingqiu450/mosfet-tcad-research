import devsim as ds, mosfet_deck as M
M.build()
I={}
for dv in [2e-4,-2e-4,5e-4,-5e-4]:
    ds.circuit_alter(name='V1',value=dv); M.solve_dc()
    I[dv]=-ds.get_circuit_node_value(solution='dcop',node='V1.I')
    print(dv,I[dv],ds.get_contact_current(device='mos',contact='body',equation='ElectronContinuityEquation'))
for h in [2e-4,5e-4]: print('G central',h,(I[h]-I[-h])/(2*h))
