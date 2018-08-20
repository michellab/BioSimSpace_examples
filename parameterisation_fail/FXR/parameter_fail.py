import BioSimSpace as BSS
import glob


def try_parametrise(filename, forcefield):
    system = BSS.IO.readMolecules(filename)
    mol = system.getMolecules()[0]
    newMol = BSS.Parameters.parameterise(mol, forcefield).getMolecule()
    if newMol is None:
        print('ERROR: System was not parametrised for file: ' + filename)
    else:
        print('SUCCESS! We have parameters for ' + filename)
        return newMol

if __name__ == '__main__':

    ligs = 'ligand.pdb'
    newMol = try_parametrise(ligs, 'gaff')