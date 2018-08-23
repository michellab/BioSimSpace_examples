import BioSimSpace as BSS


if __name__ == '__main__':

    # problematic file
    filename = 'FXR_84_BM2.mol2'

    # loading the molecule
    mol0 = BSS.IO.readMolecules(filename).getMolecules()[0]

    mapping = BSS.Align.matchAtoms(mol0, mol0)

    for idx0, idx1 in mapping.items():
        atom0 = mol0._getSireMolecule().atom(idx0)
        atom1 = mol0._getSireMolecule().atom(idx1)
        map_string = "%s --> %s" % (atom0, atom1)
        print (map_string)
