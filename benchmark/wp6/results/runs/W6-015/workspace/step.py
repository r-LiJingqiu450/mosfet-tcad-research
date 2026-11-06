import devsim as ds, mosfet_deck as M
M.build()
I=lambda: ds.get_circuit_node_value(node='V1.I', solution='dcop')
for dv in (1e-6,1e-8):
    ds.circuit_alter(name='V1', value=dv)
    out=[]
    for k in range(5):
        M._solve_one(); out.append(I()/dv)
    print('dv',dv,'I/dv per Newton iter', out)
    ds.circuit_alter(name='V1', value=0.); M.solve_dc()
for f in (1e-3,1.0,1e3):
    ds.solve(type='ac', frequency=f)
    print('AC',f, ds.get_circuit_node_value(node='V1.I', solution='ssac_real'), ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'))
