#
import BioSimSpace as BSS

m0 = BSS.IO.readMolecules('../small-molecules/methane/methane.pdb').getMolecules()[0]
m0 = BSS.Parameters.parameterise(m0, "gaff2").getMolecule()
BSS.IO.saveMolecules("m0_vac", m0, ["prm7", "rst", "pdb"])

m1 = BSS.IO.readMolecules('../small-molecules/toluene/toluene.pdb').getMolecules()[0]
m1 = BSS.Parameters.parameterise(m1, "gaff2").getMolecule()
BSS.IO.saveMolecules("m1_vac", m1, ["prm7", "rst", "pdb"])

mapping = BSS.Align.matchAtoms(m0, m1, verbose=True, scoring_function='RMSD align')#, match_light=False)

print (mapping)
m0 = BSS.Align.rmsdAlign(m0, m1, mapping)
aligned_mols = m0 + m1
BSS.IO.saveMolecules('aligned.pdb', aligned_mols, "PDB")

map0 = { "coordinates" : "coordinates0",
        "charge" : "charge0",
        "element" : "element0" }

map1 = { "coordinates" : "coordinates1",
        "charge" : "charge1",
        "element" : "element1" }

merged = BSS.Align.merge(m0, m1, mapping=mapping)
BSS.IO.saveMolecules('merged0.pdb', merged,'PDB', map=map0)
BSS.IO.saveMolecules('merged1.pdb', merged,'PDB', map=map1)
