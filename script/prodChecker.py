import os
import pybel
import rmgpy.molecule
from rmgpy.molecule.converter import from_ob_mol


def toRMGmol(OBMol):
    rmg_mol = from_ob_mol(rmgpy.molecule.molecule.Molecule(), OBMol)
    return rmg_mol


def readXYZ(path):
    mol = next(pybel.readfile('xyz', path))
    return mol.OBMol


def molChecker():
    dirs = os.listdir('/mnt/d/reactions')
    for i in dirs:
        ssm_product = os.path.join(os.path.join(
            '/mnt/d/reactions', i), 'ssm_product.xyz')
        OBMol_1 = readXYZ(ssm_product)
        rmg_mol_1 = toRMGmol(OBMol_1)
        ard_product = os.path.join(os.path.join(
            '/mnt/d/reactions', i), 'product.xyz')
        OBMol_2 = readXYZ(ard_product)
        rmg_mol_2 = toRMGmol(OBMol_2)
        if rmg_mol_1.to_inchi_key() != rmg_mol_2.to_inchi_key():
            print('not match', i)


molChecker()
