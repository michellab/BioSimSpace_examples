import os
import csv
from rdkit import Chem
from rdkit.Chem import rdmolops, Descriptors, rdDepictor
from rdkit.Chem.rdMolDescriptors import Properties
from rdkit.Chem.Draw import IPythonConsole, MolsToGridImage
import py3Dmol
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from cresset import flare
import urllib.request


def filter_extract_mol(row, headers_dict):
    relation_idx = headers_dict['RELATION']
    ro5_idx = headers_dict['NUM_RO5_VIOLATIONS']
    pchembl_idx = headers_dict['PCHEMBL_VALUE']
    try:
        ro5_violations = int(row[ro5_idx])
    except:
        return None
    try:
        float(row[pchembl_idx])
    except:
        return None
    if (row[relation_idx] != '=' or ro5_violations > 0):
        return None
    mol = Chem.MolFromSmiles(row[smiles_idx])
    if (not mol):
        return None
    mols = list(rdmolops.GetMolFrags(mol, asMols = True))
    if (not mols):
        return None
    mols.sort(reverse = True, key = lambda m: m.GetNumAtoms())
    mol = mols[0]
    molWt = Descriptors.MolWt(mol)
    if (molWt < 300 or molWt > 400):
        return None
    return mol

def view_3d(ligand):
    p3 = py3Dmol.view(width = 400, height = 300)
    p3.addModel(Chem.MolToMolBlock(ligand.to_rdmol()), 'sdf')
    p3.setStyle({'stick':{}})
    p3.zoomTo()
    return p3

if __name__ == '__main__':
    table = []
    with open('/Users/toni_brain/Projects/git/BioSimSpace_examples/Flare/bioactivity-18_11_10_51.txt', encoding='utf-8') as hnd:
        rows = csv.reader(hnd, delimiter='\t')
        for i, row in enumerate(rows):
            if (i == 0):
                headers = list(row)
            else:
                table.append(list(row))

    headers_dict = {h: i for i, h in enumerate(headers)}
    smiles_idx = headers_dict['CANONICAL_SMILES']
    pchembl_idx = headers_dict['PCHEMBL_VALUE']
    chemblid_idx = headers_dict['CMPD_CHEMBLID']

    p = flare.Project()
    if (flare.main_window()):
        flare.main_window().project = p

    smiles_set = set()
    props = Properties()
    n_vec = props.GetPropertyNames()

    rdk_mols = []

    for row in table:
        mol = filter_extract_mol(row, headers_dict)
        if (mol is None):
            continue
        row[smiles_idx] = Chem.MolToSmiles(mol)
        if (row[smiles_idx] in smiles_set):
            continue
        smiles_set.add(row[smiles_idx])
        mol.SetProp('_Name', row[chemblid_idx])
        rdmolops.AssignStereochemistry(mol)
        p_vec = props.ComputeProperties(mol)
        too_flexible = False
        for name, value in zip(n_vec, p_vec):
            if (name == 'NumRotatableBonds' and value > 5):
                too_flexible = True
                break
            mol.SetProp('RDK{0:s}'.format(name), '{0:.3f}'.format(value))
        if (too_flexible):
            continue
        for i in range(pchembl_idx + 1):
            mol.SetProp(headers[i], row[i])
        rdDepictor.Compute2DCoords(mol)
        p.ligands.append(mol)
        rdk_mols.append(mol)

    for l in p.ligands:
        l.add_hydrogens()
        l.minimize()
        filename = '/Users/toni_brain/Projects/git/BioSimSpace_examples/Flare/Thrombin_ligands/'+l.title+str('.mol2')
        l.write_file(filename)

        







