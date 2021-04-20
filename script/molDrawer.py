import os
import pybel
import rmgpy.molecule
from rmgpy.molecule.converter import from_ob_mol
from rmgpy.molecule.draw import MoleculeDrawer


def toRMGmol(OBMol):
    rmg_mol = from_ob_mol(rmgpy.molecule.molecule.Molecule(), OBMol)
    return rmg_mol


def readXYZ(path):
    mol = next(pybel.readfile('xyz', path))
    return mol.OBMol


def molDrawer():
    dirs = os.listdir('/mnt/d/reactions')
    for i in dirs:
        dir_path = os.path.join(os.path.join(
            '/mnt/d/reactions', i), 'ssm_product.xyz')
        OBMol = readXYZ(dir_path)
        rmg_mol = toRMGmol(OBMol)
        _path = '/mnt/d/molecules/{}.png'.format(i)
        MoleculeDrawer().draw(rmg_mol, file_format='png', target=_path)


def generate_ssm_product_xyz():
    opt_file = []
    reactions_path = '/mnt/d/reactions'
    reactions_path_list = os.listdir(reactions_path)
    for i in reactions_path_list:
        target_path = os.path.join(reactions_path, i)
        path_list = os.listdir(os.path.join(target_path, 'SSM'))
        path_list.sort()
        for filename in path_list:
            if filename.startswith('opt'):
                opt_file.append(filename)
        product_xyz = os.path.join(os.path.join(
            target_path, 'SSM'), opt_file[-1])
        with open(product_xyz, 'r') as f:
            lines = f.readlines()
            for i in reversed(lines):
                a = i.split()
                if len(a) == 1:
                    idx = lines.index(i)
                    break
        parent_ssm_product_path = os.path.join(os.path.abspath(
            os.path.join(os.path.join(target_path, 'SSM'), '../')), 'ssm_product.xyz')
        with open(parent_ssm_product_path, 'w') as q:
            q.write('{}\n{}'.format(lines[idx-1], ''.join(lines[idx+1:])))

# molDrawer()
# generate_ssm_product_xyz()
