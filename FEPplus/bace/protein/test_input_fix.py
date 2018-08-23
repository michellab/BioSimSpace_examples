import BioSimSpace as BSS
protein = BSS.IO.readMolecules("protein_only.pdb").getMolecules()[0]
process = BSS.Parameters.parameterise(protein, "ff14SB")
process.getOutput("param")
