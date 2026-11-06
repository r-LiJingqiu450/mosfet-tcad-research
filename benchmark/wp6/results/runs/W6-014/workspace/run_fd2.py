import sys, time, json
import numpy as np
import devsim as ds
import mosfet_deck as M

EXT = sys.argv[1] == 'ext'
if EXT:
    for p in ('extended_solver', 'extended_model', 'extended_equation'):
        ds.set_parameter(name=p, value=True)

t0 = time.time()
info = M.build()
print('build', info, time.time() - t0, flush=True)

for v in ('x', 'y'):
    ds.edge_from_node_model(device=M.DEV, region=M.SI, node_model=v)
x0 = np.array(ds.get_edge_model_values(device=M.DEV, region=M.SI, name='x@n0'))
x1 = np.array(ds.get_edge_model_values(device=M.DEV, region=M.SI, name='x@n1'))
ec = np.array(ds.get_edge_model_values(device=M.DEV, region=M.SI, name='EdgeCouple'))
cuts = [-1.5e-5, 0.1e-4 * 0.0 + 1e-6, 0.5e-4, 0.99e-4, 1.15e-4]
masks = []
for xc in cuts:
    s = np.where((x0 < xc) & (x1 > xc), 1.0, np.where((x1 < xc) & (x0 > xc), -1.0, 0.0))
    masks.append(s * ec)


def currents():
    out = {}
    for c in ('source', 'drain', 'body'):
        out[c] = ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation')
    out['V1.I'] = ds.get_circuit_node_value(solution='dcop', node='V1.I')
    j = np.array(ds.get_edge_model_values(device=M.DEV, region=M.SI, name='ElectronCurrent'))
    out['cuts'] = [float(np.sum(m * j)) for m in masks]
    return out


res = {}
res[0.0] = currents()
print(0.0, res[0.0], flush=True)

for f in (1e-6, 1e-3):
    ds.solve(type='ac', frequency=f)
    re = ds.get_circuit_node_value(solution='ssac_real', node='V1.I')
    im = ds.get_circuit_node_value(solution='ssac_imag', node='V1.I')
    print('AC', f, re, im, flush=True)
    res[f'AC{f}'] = [re, im]

for dv in [float(a) for a in sys.argv[2:]]:
    ds.circuit_alter(name='V1', value=dv)
    h = M.solve_dc()
    res[dv] = currents()
    print(dv, h['iterations'], h['history'][-3:], res[dv], flush=True)
    ds.circuit_alter(name='V1', value=0.0)
    M.solve_dc()
res['zero_again'] = currents()
print('zero again', res['zero_again'])
json.dump({str(k): v for k, v in res.items()}, open(f'fd2_{sys.argv[1]}.json', 'w'), indent=1)
