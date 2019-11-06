import BioSimSpace as BSS

#load the two molecules
ethane = BSS.IO.readMolecules(BSS.IO.glob("molecules/ethane.*"))[0]
methanol = BSS.IO.readMolecules(BSS.IO.glob("molecules/methanol.*"))[0]

# now work out the maximum common substructure and align them
mapping = BSS.Align.matchAtoms(ethane, methanol)
ethane = BSS.Align.rmsdAlign(ethane, methanol, mapping)

# now we create a merged molecule that will be able to deal with intermediates
merged = BSS.Align.merge(ethane, methanol, mapping)

# Next solvate the system
solvated = BSS.Solvent.tip3p(molecule=merged, box=3*[5*BSS.Units.Length.nanometer])

# Give it the default free energy protocol
protocol = BSS.Protocol.FreeEnergy(num_lam=5, runtime=100*BSS.Units.Time.picosecond)

# Run the free energy simulation
freenrg_somd = BSS.FreeEnergy.Solvation(solvated, protocol, engine="SOMD", work_dir="freenrg_somd")

freenrg_somd.run()
freenrg_somd.analyse()
