import BioSimSpace as BSS

if __name__ == '__main__':

    # Load the molecule.
    mol0 = BSS.IO.readMolecules("FXR_84_BM2.mol2").getMolecules()[0]
    mol1 = BSS.IO.readMolecules("FXR_84_BM2_shifted.mol2").getMolecules()[0]

    # Get the best mapping.
    mapping0 = BSS.Align.matchAtoms(mol0, mol0)
    mapping1 = BSS.Align.matchAtoms(mol1, mol1)

    # Make sure all of the atoms have been mapped.
    print("Unshifted ==> Mapped %d out of %d atoms." % (len(mapping0), mol0.nAtoms()))
    print("Shifted   ==> Mapped %d out of %d atoms." % (len(mapping1), mol1.nAtoms()))
