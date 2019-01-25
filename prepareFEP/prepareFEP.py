import BioSimSpace as BSS
import os
#
#
#
# Version 1)
# Reads a pair of ligands A and B, solvated by water or water+protein+ions 
# Writes SOMD inputs for a perturbation between ligands A to B
#
# input requirements.
# prm7/rst7 files for A and for B
# output
# SOMD prm7/rst7 + pert files for AtoB
#
# TODO
# - convert in notebook
# - log mapping details
# - interface to parse input
# - add error checking
# - make mappings consistent between forwwards and backwards options
# - test Gromacs support
# - give option to specify output for desired MD package


#inputA_top = "input/1a.parm7"
#inputA_rst = "input/1a_equilibrated.rst7"
#inputB_top = "input/5.parm7"
#inputB_rst = "input/5_equilibrated.rst7"
inputA_top = "input/THROM-1a.parm7"
inputA_rst = "input/THROM-1a_equilibrated.rst7"
inputB_top = "input/THROM-5.parm7"
inputB_rst = "input/THROM-5_equilibrated.rst7"
# Load system A
systemA = BSS.IO.readMolecules([inputA_top, inputA_rst])
# Extract ligand (assumed first molecule) and unperturbed molecules
moleculesA = systemA.getMolecules()
ligA = moleculesA[0]
remainderA = moleculesA[1:]
# Load system B
systemB = BSS.IO.readMolecules([inputB_top, inputB_rst])
moleculesB = systemB.getMolecules()
ligB = moleculesB[0]
remainderB = moleculesB[1:]
# Perform mapping between ligA and ligB
mapping = BSS.Align.matchAtoms(ligA, ligB)
#
# @ Optimisation. Retain box with largest number of atoms for consistency between A->B and B->A mappings
#
# Align ligB to ligA based on the mapping.
ligB = BSS.Align.rmsdAlign(ligB, ligA, mapping)
# Merge the two ligands based on the mapping.
merged = BSS.Align.merge(ligA, ligB, mapping)
# Create a composite system
systemA.removeMolecules(ligA)
systemA.addMolecules(merged)
#system = merged + remainderA
protocol = BSS.Protocol.FreeEnergy(runtime = 2*BSS.Units.Time.femtosecond, num_lam=3)
# Generate package specific input
process = BSS.Process.Somd(systemA, protocol)
process.getOutput()
cmd = "unzip -o somd.zip"
os.system(cmd)
pert = inputA_top.split(".")[0] + "_to_" + inputB_top.split(".")[0]
cmd = "mv somd.pert output/%s.pert ; mv somd.prm7 output/%s.prm7 ; mv somd.rst7 output/%s.rst7 ; rm somd.zip ; rm somd.cfg ; rm somd.err; rm somd.out" % (pert,pert,pert)
print (cmd)
os.system(cmd)
