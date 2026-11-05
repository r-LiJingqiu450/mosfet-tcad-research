import sys, time, json
import devsim as ds
import mosfet_deck as M

t0 = time.time()
info = M.build()
print('build', info, time.time() - t0, flush=True)


def currents():
    out = {}
    for c in ('source', 'drain', 'body'):
        out[c] = ds.get_contact_current(device=M.DEV, contact=c, equation='ElectronContinuityEquation')
    out['V1.I'] = ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return out


res = {}
res[0.0] = currents()
print(0.0, res[0.0], flush=True)

# small-signal AC at very low frequency
for f in (1e-3, 1.0):
    ds.solve(type='ac', frequency=f)
    re = ds.get_circuit_node_value(solution='ssac_real', node='V1.I')
    im = ds.get_circuit_node_value(solution='ssac_imag', node='V1.I')
    print('AC', f, re, im, flush=True)

for dv in [float(a) for a in sys.argv[1:]]:
    ds.circuit_alter(name='V1', value=dv)
    h = M.solve_dc()
    res[dv] = currents()
    print(dv, h['iterations'], res[dv], flush=True)
    ds.circuit_alter(name='V1', value=0.0)
    M.solve_dc()
json.dump({str(k): v for k, v in res.items()}, open('fd_results.json', 'w'), indent=1)
