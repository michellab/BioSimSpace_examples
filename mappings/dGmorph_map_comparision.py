import BioSimSpace as BSS

from os.path import basename
import itertools
import sys


def parse_dGmorph(filename):

    print('Parsing the morph file....')

    fh = open(filename, 'r')
    lines = fh.readlines()
    generating_morphs = False
    morph_pairs = {}
    perturbation = None
    mapping_list = []

    ## looping over lines in file to extract relevant mapping
    for line in lines:
        if generating_morphs:
            curr_perturbation = perturbation
            if line.startswith('Entering'):
                if '_perturbations' in line:
                    perturbation = (line.split('/')[-1].strip())
                    if curr_perturbation is None:
                        curr_perturbation = perturbation
            if line.startswith('AtomName'):
                l = line.strip()
                temp = l.split('\'')
                mapping_list.append([temp[1],temp[3]])
            if curr_perturbation != perturbation:
                morph_pairs[curr_perturbation] = mapping_list
                mapping_list = []

        #lines after this will contain The atom mapping
        if line.startswith('Morphs will be generated'):
            generating_morphs = True
    morph_pairs[perturbation] = mapping_list

    print('Parsing of the morph file done...')

    return morph_pairs



def test_mapping(file0, file1, count, output, verbose=False, timeout = '5 secs'):

    BSS_mapping = []

    if verbose:
        print("Loading the molecules...")

    #print("%05d : %s --> %s" % (count, basename(file0), basename(file1)))

    # Load each file and grab the first (only) molecule.
    mol0 = BSS.IO.readMolecules(file0).getMolecules()[0]
    mol1 = BSS.IO.readMolecules(file1).getMolecules()[0]

    if verbose:
        print("Performing the mapping...")

    # Find the best MCS mapping.
    mapping = BSS.Align.matchAtoms(mol0, mol1, verbose=verbose, timeout=BSS.Types.Time(timeout))

    # Write the mapping to file and list
    with open(output + ".mapping", "w") as file:
        for idx0, idx1 in mapping.items():
            atom0 = mol0._getSireMolecule().atom(idx0)
            atom1 = mol1._getSireMolecule().atom(idx1)
            map_string = "%s --> %s\n" % (atom0, atom1)
            temp = map_string.split(' ')
            BSS_mapping.append([temp[1],temp[7]])
            file.write("%s --> %s\n" % (atom0, atom1))

    # Align mol0 to mol1.
    mol0 = BSS.Align.rmsdAlign(mol0, mol1, mapping)

    # Create a combined system from the aligned molecules.
    system = mol0 + mol1

    # Write the system to file.
    BSS.IO.saveMolecules(output, system, "PDB")
    return BSS_mapping


def compare_mapping(FESetup_mapping, BioSimSpace_mapping, ligpair, verbose=False):

    if len(FESetup_mapping) != len(BioSimSpace_mapping):
        print('Mappings are not of the same length for ligand pair %s' %ligpair)
        print('FESetup length is %d' % len(FESetup_mapping))
        print('BSS length is %d' % len(BioSimSpace_mapping))
        if verbose:
            for i in range(len(FESetup_mapping)):
                if i < len(BioSimSpace_mapping):
                    print(' Fesetup: atom1 %s --> atom2 %s  | BSS: atom1 %s --> atom2 %s' %(FESetup_mapping[i][0], FESetup_mapping[i][1], BioSimSpace_mapping[i][0], BioSimSpace_mapping[i][1]))
                else:
                    print(' Fesetup: atom1 %s --> atom2 %s ' %(FESetup_mapping[i][0], FESetup_mapping[i][1]))
        return False
    print(FESetup_mapping)
    print(BioSimSpace_mapping)
    return FESetup_mapping is BioSimSpace_mapping


def write_mapping_FESetup(output, atoms):
    with open(output + ".mapping", "w") as file:
        for a in atoms:
            file.write("%s --> %s\n" % (a[0], a[1]))


if __name__ == '__main__':

    # Getting info from the commandline
    if len(sys.argv) < 3:
        print('Usage: dGmorph_map_comparision.py dGmorph.log /dir/to/_ligands/ verbose=False')

    dGmorph = sys.argv[1]
    inputdirectory = sys.argv[2]
    verbose = False

    if len(sys.argv) >3:
        verbose = bool(sys.argv[3])
        print('Verbosity set to True')

    FESetup_mapping = parse_dGmorph(dGmorph)

    # now generate the mappings for BSS
    BSS_mapping = {}

    # counting the ligand pairs
    count = 0

    # looping over all keys in the FESetup mapping to generate the same maps


    for lig_pair, lig_pair_map in FESetup_mapping.items():
        pair = lig_pair.split('~')
        file0 = inputdirectory + pair[0]+'/gaff.mol2'
        file1 = inputdirectory+ pair[1]+'/gaff.mol2'
        output_BSS = "BSS_%02d" % count + "_" + pair[0] + "_" + pair[1]
        mapping = test_mapping(file0, file1, count, output_BSS, verbose)
        BSS_mapping[lig_pair] = mapping


        #comparing the two mappings
        if compare_mapping(FESetup_mapping[lig_pair], mapping, lig_pair, verbose):
            print('Mapping for %s perturbation is identical' %lig_pair)
        else:
            print('There is a problem with the mapping of %s ' %lig_pair)
            print(' Look at the mapping files to see what the problem may be!')
            output_FES = "FES_%02d" % count + "_" + pair[0] + "_" + pair[1]
            write_mapping_FESetup(output_FES, FESetup_mapping[lig_pair])

        count = count + 1