import mosfet_deck as M, devsim as ds, json
info = M.build(); print(info)
def cur():
    return {c: ds.get_contact_current(device='mos', contact=c, equation='ElectronContinuityEquation') for c in ('source','drain','body')}, ds.get_circuit_node_value(solution='dcop', node='V1.I')
print(0, cur())
res={}
for dv in (1e-4,-1e-4,1e-3,-1e-3):
    ds.circuit_alter(name='V1', value=dv); M.solve_dc()
    res[dv]=cur(); print(dv, res[dv])
ds.circuit_alter(name='V1', value=0.0); M.solve_dc()
for d in (1e-4,1e-3):
    print('central', d, (res[d][1]-res[-d][1])/(2*d), (res[d][0]['drain']-res[-d][0]['drain'])/(2*d))
ds.solve(type='ac', frequency=1e-3)
print('ac', ds.get_circuit_node_value(solution='ssac_real', node='V1.I'), ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
