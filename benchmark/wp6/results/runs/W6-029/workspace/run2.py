import mosfet_deck as M, devsim as ds
M.build()
for n in ['PotentialNodeCharge:Potential','PotentialNodeCharge:Electrons']:
    print(n, ds.get_node_model_list(device='mos',region='si').__contains__(n))
print(ds.get_equation_command(device='mos', region='si', name='PotentialEquation') if hasattr(ds,'get_equation_command') else '')
for v in [1e-4]:
    ds.circuit_alter(name='V1', value=v)
    for k in range(8):
        info=ds.solve(type='dc', absolute_error=1e300, relative_error=1e300, maximum_iterations=1, solver_type='direct', info=True)
        e=M._errors(info); print(k, {kk:(vv if not isinstance(vv,dict) else vv['abs']) for kk,vv in e.items()}, ds.get_circuit_node_value(solution='dcop', node='V1.I'))
