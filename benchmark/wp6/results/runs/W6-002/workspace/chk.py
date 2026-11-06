import devsim as ds, mosfet_deck as M
M.build()
for r in (M.SI,):
    for nm in ('ElectronCurrent:Electrons@n0','ElectronCurrent:Electrons@n1','ElectronCurrent:Potential@n0','ElectronCurrent:Potential@n1','PotentialNodeCharge:Potential','PotentialNodeCharge:Electrons'):
        try: print(nm, ds.get_edge_model_values if 0 else '', ds.get_equation_command if 0 else '')
        except Exception as e: print(e)
ds.circuit_alter(name='V1', value=1e-3)
h=M.solve_dc()
for x in h['history']: print(x)
