import devsim as ds, mosfet_deck as M
M.build_mesh(); M.build_physics()
for n in ('ElectronCurrent:Potential@n0','ElectronCurrent:Electrons@n0','ElectronCurrent:Electrons@n1'):
    print(n, '=', ds.edge_model(device=M.DEV, region=M.SI, name=n+'_probe', equation=f'simplify(diff(ElectronCharge*mu_n*EdgeInverseLength*V_t*kahan3(Electrons@n1*Bern01,  Electrons@n1*vdiff,  -Electrons@n0*Bern01),{n.split(":")[1]}))'))
print(ds.node_model(device=M.DEV, region=M.SI, name='pp', equation='simplify(diff(-ElectronCharge*kahan3(Holes, -Electrons, NetDoping),Potential))'))
