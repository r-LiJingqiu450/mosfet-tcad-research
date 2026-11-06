import devsim as ds, mosfet_deck as M
M.build()
def cur():
    return [ds.get_circuit_node_value(solution='dcop', node='V1.I')]+[ds.get_contact_current(device='mos', contact=c, equation='ElectronContinuityEquation') for c in ('drain','source','body')]
print('R0', cur())
res={}
for dv in [1e-2,-1e-2,3e-3,-3e-3]:
    ds.circuit_alter(name='V1', value=dv); M.solve_dc(); res[dv]=cur(); print('R',dv,res[dv], sum(res[dv][1:]))
for dv in [1e-2,3e-3]:
    print('G', dv, [(a-b)/(2*dv) for a,b in zip(res[dv],res[-dv])])
