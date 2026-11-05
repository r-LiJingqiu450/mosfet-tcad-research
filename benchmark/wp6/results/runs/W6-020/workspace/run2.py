import mosfet_deck as M, devsim as ds, json, numpy as np
info = M.build(); print(info)
I=lambda: ds.get_circuit_node_value(solution='dcop', node='V1.I')
print('I0', I())
for dv in (1e-5,-1e-5, 1e-2, -1e-2):
    ds.circuit_alter(name='V1', value=dv); r=M.solve_dc()
    print(dv, I(), I()/dv, [(f"{h['dpsi']:.1e}",f"{h['dphin']:.1e}") for h in r['history']])
    for c in ('source','drain','body'):
        print('  ',c, ds.get_contact_current(device='mos', contact=c, equation='ElectronContinuityEquation'), ds.get_contact_charge(device='mos', contact=c, equation='PotentialEquation'))
