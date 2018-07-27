import BioSimSpace as BSS

from os.path import basename

import itertools
import sys

def test_mapping(file0, file1, count, output, verbose=False):

    if verbose:
        print("Loading the molecules...")

    print("%05d : %s --> %s" % (count, basename(file0), basename(file1)))

    # Load each file and grab the first (only) molecule.
    mol0 = BSS.IO.readMolecules(file0).getMolecules()[0]
    mol1 = BSS.IO.readMolecules(file1).getMolecules()[0]

    if verbose:
        print("Performing the mapping...")

    # Find the best MCS mapping.
    mapping = BSS.Align.matchAtoms(mol0, mol1, verbose=verbose)

    # Write the mapping to file.
    with open(output + ".mapping", "w") as file:
        for idx0, idx1 in mapping.items():
            atom0 = mol0._getSireMolecule().atom(idx0)
            atom1 = mol1._getSireMolecule().atom(idx1)

            file.write("%s --> %s\n" % (atom0, atom1))

    # Align mol0 to mol1.
    mol0 = BSS.Align.rmsdAlign(mol0, mol1, mapping)

    # Create a combined system from the aligned molecules.
    system = mol0 + mol1

    # Write the system to file.
    BSS.IO.saveMolecules(output, system, "PDB")

if __name__ == "__main__":

    if len(sys.argv) == 2:
        ligand_dir = sys.argv[1]
    else:
        ligand_dir = "."

    ligand_list = BSS.IO.glob("%s/*.mol2" % ligand_dir)

    count = 0

    # Loop over all pairs of ligands.
    for pair in itertools.combinations(ligand_list, 2):
        file0 = pair[0]
        file1 = pair[1]
        name_0 = basename(file0).split(".")[0]
        name_1 = basename(file1).split(".")[0]
        output = "%05d" % count + "_" + name_0 + "_" + name_1

        try:
            test_mapping(file0, file1, count, output, False)
        except:
            print ("ERROR: mapping problem: %s, %s" % (file0, file1))

        count = count +1
