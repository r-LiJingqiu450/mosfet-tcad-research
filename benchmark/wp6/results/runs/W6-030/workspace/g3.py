import devsim as ds, mosfet_deck as M
M.build()
for dv in [1e-3,1e-4,1e-5]:
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc()
    print('DV',dv, r['history'][:3])
    ds.circuit_alter(name='V1', value=0.); M.solve_dc()
