import mosfet_deck as M, devsim as ds
print(M.build())
ds.circuit_alter(name='V1', value=1e-3)
r=M.solve_dc()
for h in r['history']: print(h)
