import devsim as ds
import mosfet_deck as M
M.build()
for f in [0.0, 1.0, 1e3]:
    try:
        ds.solve(type='ac', frequency=f)
        print('f',f,'ReI',ds.get_circuit_node_value(solution='ssac_real', node='V1.I'),
              'ImI',ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'),flush=True)
    except Exception as e:
        print('f',f,'err',e)
print('I0', ds.get_circuit_node_value(solution='dcop', node='V1.I'),
      ds.get_contact_current(device='mos',contact='drain',equation='ElectronContinuityEquation'),
      ds.get_contact_current(device='mos',contact='source',equation='ElectronContinuityEquation'))
