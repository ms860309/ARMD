# standard library imports
import numpy as np
from openbabel import pybel as pb
from openbabel import openbabel as ob
from subprocess import Popen, PIPE

import shutil
import sys
import os
from os import path
sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'code/ard'))
sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'code/mol_graph'))
sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), 'database'))

# local application imports
from gen3D import Molecule
import gen3D

def readXYZ(xyz, bonds=None):
    # extract molecule information from xyz
    mol = next(pb.readfile('xyz', xyz))
    # Manually give bond information
    # (Because in metal system the bond information detect by openbabel usually have some problem)
    m = Molecule(pb.ob.OBMol())
    obmol = m.OBMol
    obmol.BeginModify()
    for atom in mol:
        coords = [coord for coord in atom.coords]
        atomno = atom.atomicnum
        obatom = ob.OBAtom()
        obatom.thisown = 0
        obatom.SetAtomicNum(atomno)
        obatom.SetVector(*coords)
        obmol.AddAtom(obatom)
        del obatom

    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondOrder())
             for bond in pb.ob.OBMolBondIter(mol.OBMol)]
    bonds.extend([(12,14,1), (12,15,1), (12,16,1), (12,17,1), (12,13,1), (17,23,1),(16,23,1)])

    for bond in bonds:
        obmol.AddBond(bond[0], bond[1], bond[2])

    # obmol.PerceiveBondOrders()
    # obmol.SetTotalCharge(int(mol.charge))
    # obmol.Center()
    obmol.EndModify()
    mol_obj = gen3D.Molecule(obmol)
    return mol_obj


def extract_bonds(bonds):
    with open(bonds, 'r') as f:
        lines = f.read()
    lines = eval(lines)
    return lines


def extract_fixed_atom_index(fixed_atoms):
    with open(fixed_atoms, 'r') as f:
        lines = f.read()
    lines = eval(lines)
    return lines


def getE(tmpdir, target='reactant'):
    """
    Here the energy is Eh (hartree)
    """
    input_path = path.join(tmpdir, f'{target}.xyz')
    with open(input_path, 'r') as f:
        lines = f.readlines()
    HeatofFormation = lines[1].strip().split()[1]
    return HeatofFormation


def check_bond_length(product, add_bonds):
    """
    Use reactant coordinate to check if the add bonds's bond length is too long.
    Return a 'list of distance'.
    """
    coords = [atom.coords for atom in product]
    atoms = tuple(atom.atomicnum for atom in product)
    coords = [np.array(coords).reshape(len(atoms), 3)]

    dist = []
    for bond in add_bonds:
        coord_vect_1 = coords[0][bond[0]]
        coord_vect_2 = coords[0][bond[1]]
        diff = coord_vect_1 - coord_vect_2
        dist.append(np.linalg.norm(diff))

    if dist == []:
        dist = [0]
    return float(max(dist))

def xtb_get_H298(reactant_mol, product_mol, fixed_atoms, reactant_path):
    """
    Create a directory folder called "tmp" for mopac calculation
    Create a input file called "input.mop" for mopac calculation
    """

    tmpdir = path.join(reactant_path, 'tmp')
    reactant_path = path.join(tmpdir, 'reactant.xyz')
    product_path = path.join(tmpdir, 'product.xyz')

    reac_geo, prod_geo = gen_geometry(reactant_mol, product_mol, fixed_atoms)

    if path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.mkdir(tmpdir)
    os.chdir(tmpdir)
    reac_geo = f"{len(str(reac_geo).splitlines())}\n\n{reac_geo}"
    with open(reactant_path, 'w') as f:
        f.write(reac_geo)

    runXTB(tmpdir, 'reactant')
    reactant_energy = getE(tmpdir, 'reactant')

    prod_geo = f"{len(str(prod_geo).splitlines())}\n\n{prod_geo}"
    with open(product_path, 'w') as f:
        f.write(prod_geo)

    runXTB(tmpdir, 'product')
    product_energy = getE(tmpdir, 'product')

    return float(reactant_energy), float(product_energy)


def runXTB(tmpdir, target='reactant', constraint=True):
    outname = f'{target}.xyz'
    input_path = path.join(tmpdir, outname)

    output_path = path.join(tmpdir, 'xtbopt.xyz')
    config_path = path.join(path.dirname(path.dirname(tmpdir)), 'config')
    constraint_path = path.join(config_path, 'xtb_constraint.inp')
    if not path.exists(constraint_path):
        config_path = path.join(path.dirname(path.dirname(path.dirname(tmpdir))), 'config')
        constraint_path = path.join(config_path, 'xtb_constraint.inp')
    new_output_path = path.join(tmpdir, outname)
    if constraint == None:
        p = Popen(['xtb', input_path, '--gfn', '2', '--opt', 'tight'])
        p.wait()
        shutil.move(output_path, new_output_path)
    else:
        p = Popen(['xtb', '--opt', 'tight', '--gfn', '2', '--input', constraint_path, input_path])
        p.wait()
        shutil.move(output_path, new_output_path)



def gen_geometry(reactant_mol, product_mol, fixed_atoms):

    # Set up constraint
    constraint_forcefield = ob.OBFFConstraints()
    constraint_forcefield.AddAtomConstraint(1)
    ff = pb.ob.OBForceField.FindForceField('uff')
    Hatom = gen3D.readstring('smi', '[H]')
    

    # reactant_bonds = [tuple(sorted((bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1)) + [bond.GetBondOrder()])
    #                     for bond in pb.ob.OBMolBondIter(reactant_mol.OBMol)]
    # reactant_bonds = tuple(sorted(reactant_bonds))
    # print(reactant_bonds)

    reactant_mol.gen3D(fixed_atoms, forcefield='uff',
                    method='ConjugateGradients', make3D=False)

    product_mol.gen3D(fixed_atoms, forcefield='uff',
                    method='ConjugateGradients', make3D=False)

    #print(product_mol.write('inchikey'))
    arrange3D = gen3D.Arrange3D(
        reactant_mol, product_mol, fixed_atoms, cluster_bond = [(12,14,1), (12,15,1), (12,16,1), (12,17,1), (12,13,1), (17,23,1),(16,23,1)])
    msg = arrange3D.arrangeIn3D()
    if msg != '':
        print(msg)

    ff.Setup(Hatom.OBMol, constraint_forcefield)
    ff.SetConstraints(constraint_forcefield)
    reactant_mol.gen3D(fixed_atoms, forcefield='uff', method='ConjugateGradients', make3D=False)
    ff.Setup(Hatom.OBMol, constraint_forcefield)
    ff.SetConstraints(constraint_forcefield)
    product_mol.gen3D(fixed_atoms, forcefield='uff', method='ConjugateGradients', make3D=False)
    ff.Setup(Hatom.OBMol, constraint_forcefield)  # Ensures that new coordinates are generated for next molecule (see above)
    ff.SetConstraints(constraint_forcefield)

    product_geometry = product_mol.toNode()
    reactant_geometry = reactant_mol.toNode()
    return reactant_geometry, product_geometry


xyz_path = './reactant.xyz'
fixed_atom_path = './fixed_atoms.txt'

fixed_atoms = extract_fixed_atom_index(fixed_atom_path)
reactant = readXYZ(xyz_path)
atoms = tuple(atom.atomicnum for atom in reactant)
rb = [(22, 27, 1), (20, 26, 1), (16, 18, 1), (22, 29, 1), (20, 28, 1), (16, 22, 1), (22, 34, 1), (15, 22, 1), (20, 30, 1), (14, 20, 1), 
(0, 7, 1), (21, 32, 1), (21, 33, 1), (0, 6, 1), (0, 1, 1), (0, 5, 1), (11, 12, 1), (12, 21, 1), (11, 14, 1), (15, 17, 1), (3, 8, 1), (1, 2, 2), 
(11, 13, 1), (21, 31, 1), (1, 3, 1), (3, 9, 1), (3, 4, 1), (4, 10, 1), (13, 19, 1), (19, 24, 1), (19, 23, 1), (19, 25, 1), (11, 15, 1), (11, 16, 1)]


# Isomerization Hydroxyacetone --> 2-Hydroxypropanal (Acetol)
product_bonds = [(22, 27, 1), (20, 26, 1), (16, 18, 1), (22, 29, 1), (20, 28, 1), (16, 22, 1), (22, 34, 1), (15, 22, 1), (20, 30, 1), (14, 20, 1), 
                (0, 7, 1), (21, 32, 1), (21, 33, 1), (0, 6, 1), (0, 1, 1), (0, 5, 1), (11, 12, 1), (12, 21, 1), (11, 14, 1), (15, 17, 1), (3, 8, 1), (1, 2, 2), 
                (11, 13, 1), (21, 31, 1), (1, 3, 1), (3, 9, 1), (3, 4, 1), (13, 19, 1), (19, 24, 1), (19, 23, 1), (19, 25, 1), (11, 15, 1), (11, 16, 1), (4, 11, 1), (10, 14, 1)]


"""

# Fail
product_bonds = ((0, 1, 1), (0, 5, 1), (0, 6, 1), (0, 7, 1), (1, 2, 2), (1, 3, 1), (2, 11, 1), (3, 4, 1), (3, 8, 1), (3, 9, 1), (4, 18, 1), (10, 13, 1), (11, 12, 1), (11, 13, 1), (11, 14, 1), (11, 15, 1), (11, 16, 1), (12, 21, 1), (13, 19, 1), (14, 20, 1), (15, 17, 1), (15, 22, 1), (16, 22, 1), (19, 23, 1), (19, 24, 1), (19, 25, 1), (20, 26, 1), (20, 28, 1), (20, 30, 1), (21, 31, 1), (21, 32, 1), (21, 33, 1), (22, 27, 1), (22, 29, 1), (22, 34, 1))
"""
# Acrolein
"""
product_bonds = ((0, 1, 2), (0, 5, 1), (0, 6, 1), (15, 7, 1), (1, 2, 1), (1, 3, 1), (3, 4, 2), (1, 8, 1), (3, 9, 1), (2, 10, 1), (11, 12, 1), 
    (11, 13, 1), (11, 14, 1), (11, 15, 1), (11, 16, 1), (12, 21, 1), (13, 19, 1), (14, 20, 1), (2, 17, 1), (15, 22, 1), (16, 18, 1), (16, 22, 1), 
    (19, 23, 1), (19, 24, 1), (19, 25, 1), (20, 26, 1), (20, 28, 1), (20, 30, 1), (21, 31, 1), (21, 32, 1), (21, 33, 1), (22, 27, 1), (22, 29, 1), (22, 34, 1))
"""

# Prop-2-ene-1,2-diol
"""
product_bonds = ((0, 1, 2), (0, 5, 1), (0, 6, 1), (14, 7, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1), (3, 8, 1), (3, 9, 1), (4, 10, 1), (11, 12, 1),
                 (11, 13, 1), (11, 14, 1), (11, 15, 1), (11, 16, 1), (12, 21, 1), (13,
                                                                                   19, 1), (14, 20, 1), (2, 17, 1), (15, 22, 1), (16, 18, 1), (16, 22, 1),
                 (19, 23, 1), (19, 24, 1), (19, 25, 1), (20, 26, 1), (20, 28, 1), (20, 30, 1), (21, 31, 1), (21, 32, 1), (21, 33, 1), (22, 27, 1), (22, 29, 1), (22, 34, 1))
"""

product = gen3D.makeMolFromAtomsAndBonds(atoms, product_bonds, spin=reactant.spin)
product.setCoordsFromMol(reactant)

a, b = xtb_get_H298(reactant, product, fixed_atoms, os.getcwd())
print(a)
print(b)
print((b-a)*627.5095)