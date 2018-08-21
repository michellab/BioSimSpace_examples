import BioSimSpace as BSS
# Read ligand and save solvated description 
ligand = BSS.IO.readMolecules("../small-molecules/fxr_79_pdb/fxr_79.pdb").getMolecules()[0]
ligand = BSS.Parameters.parameterise(ligand, "gaff2").getMolecule()
system = BSS.Solvent.solvate("tip3p", molecule=ligand,
                                      box=3*[2.5*BSS.Units.Length.nanometer])
BSS.IO.saveMolecules("solvated_ligand", system, ["grotop", "gro87"])

# Load protein
protein = BSS.IO.readMolecules("../protein/FXR/fxr-allh.pdb").getMolecules()[0]
protein = BSS.Parameters.parameterise(protein, "ff14SB").getMolecule()

system = BSS.Solvent.solvate("tip3p", molecule=protein+ligand,
                                      box=3*[7*BSS.Units.Length.nanometer])
BSS.IO.saveMolecules("solvated_complex", system, ["prm7", "rst", "pdb"])
