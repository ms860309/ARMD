#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contains functions and classes for generating 3D geometries using Open Babel.
Also contains functionality for estimating thermo using group additivity and
RMG database values.
"""

# standard library imports
import os
from os import path
from operator import itemgetter

# third party
import numpy as np
from scipy.optimize import minimize, Bounds
from scipy.spatial import distance_matrix
import scipy.constants as spc
from openbabel import pybel
from openbabel import openbabel as ob
try:
    from rmgpy import settings
    from rmgpy.species import Species
    import rmgpy.molecule
    from rmgpy.data.thermo import ThermoDatabase
    from rmgpy.molecule.adjlist import to_adjacency_list
    from rmgpy.molecule.converter import from_ob_mol
except:
    pass
from itertools import combinations

# local application imports
from _constants import BOHR2ANG, AU2KJPERMOL, AU2KCALMOL
import _constants
import props
import util
import node
from quantum import QuantumError

ELEMENT_TABLE = props.ElementData()

###############################################################################


def readstring(format, string):
    """
    Read in a molecule from a string and convert to a :class:`Molecule` object.
    """
    mol = pybel.readstring(format, string)
    return Molecule(mol.OBMol)


def make3DandOpt(mol, forcefield='uff', make3D=True):
    """
    Generate 3D coordinates and optimize them using a force field.
    """
    if make3D:
        mol.make3D(forcefield=forcefield)
    mol.localopt(forcefield=forcefield)


def constraint_force_field(mol, freeze_index, forcefield='uff', method='ConjugateGradients', steps=5000):
    """
    An openbabel constraint force field.
    mol is OBMol object
    freeze_index is the list of atom you want to freeze (Start from 0)
    """
    # Set up constraint
    constraint_forcefield = ob.OBFFConstraints()
    for constraint_atom_index in freeze_index:
        atom_index = constraint_atom_index + 1
        constraint_forcefield.AddAtomConstraint(atom_index)

    # Set up forcefield
    ff = ob.OBForceField.FindForceField(forcefield)
    ff.Setup(mol.OBMol, constraint_forcefield)

    ff.SetConstraints(constraint_forcefield)
    if method == 'SteepestDescent':
        ff.SteepestDescent(steps)
    elif method == 'ConjugateGradients':
        ff.ConjugateGradients(steps)
    ff.GetCoordinates(mol.OBMol)

def get_ob_forces(atoms, coords, forcefield = 'mmff94', cluster_bond = [], freeze_index = []):
    xyz = prepare_xyz_string(atoms, coords)
    mol = pybel.readstring("xyz", xyz)
    m = Molecule(pybel.ob.OBMol())
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
                for bond in pybel.ob.OBMolBondIter(mol.OBMol)]
    bonds.extend(cluster_bond)
    for bond in bonds:
        obmol.AddBond(bond[0], bond[1], bond[2])
    obmol.EndModify()
    mol_obj = Molecule(obmol)

    ff = ob.OBForceField.FindForceField(forcefield)
    if freeze_index != []:
    # Set up constraint
        constraint_forcefield = ob.OBFFConstraints()
        for constraint_atom_index in freeze_index:
            atom_index = constraint_atom_index + 1
            constraint_forcefield.AddAtomConstraint(atom_index)
        ff.Setup(mol_obj.OBMol, constraint_forcefield)
        ff.SetConstraints(constraint_forcefield)
    else:
        ff.Setup(mol_obj.OBMol)

    conv_dict = {
        "kj/mol": AU2KJPERMOL,
        "kcal/mol": AU2KCALMOL,
    }
    unit = ff.GetUnit()
    #an_grad = ff.HasAnalyticalGradients()
    conv_fac = conv_dict[unit.lower()]
    energy = ff.Energy(True) / conv_fac
    grad = [ff.GetGradient(atom.OBAtom) for atom in mol_obj.atoms]
    grad = np.array([(g.GetX(), g.GetY(), g.GetZ()) for idx, g in enumerate(grad) if idx not in freeze_index])
    forces = grad.flatten() * BOHR2ANG / conv_fac
    norm_forces = np.linalg.norm(forces)
    #print(mol_obj.write('can'))
    #print(energy)
    #print(norm_forces)
    return norm_forces

def prepare_xyz_string(atoms, coords):
    coords = coords.reshape(-1, 3) * spc.value("Bohr radius") * 1e10
    coords = "\n".join(
            ["{} {:10.08f} {:10.08f} {:10.08f}".format(a, *c) for a, c in zip(atoms, coords)]
    )
    return f"{len(atoms)}\n\n{coords}"

def makeMolFromAtomsAndBonds(atoms, bonds, spin=None):
    """
    Create a new Molecule object from a sequence of atoms and bonds.
    """
    mol = Molecule(pybel.ob.OBMol())
    OBMol = mol.OBMol

    for atomicnum in atoms:
        a = pybel.ob.OBAtom()
        a.SetAtomicNum(atomicnum)
        OBMol.AddAtom(a)
    for bond in bonds:
        if len(bond) != 3:
            raise Exception(
                'Bond must be specified by two indices and a bond order')
        OBMol.AddBond(bond[0] + 1, bond[1] + 1, bond[2])

    mol.assignSpinMultiplicity()
    if spin is not None:
        OBMol.SetTotalSpinMultiplicity(spin)
    # OBMol.SetHydrogensAdded()

    return mol

def det_connectivity(total_bond, atom_index, pass_idx):
    """
    For a given atom index return its connectivity.
    The index you want to pass.
    For example, ethane (1H-2H-3H)-4C-5C-(6H-7H-8H)
    With atom_index == 4, pass index == 5, return 1,2,3
    Here index should start from 0.
    """
    atom_set = bond_search(total_bond, {atom_index}, {atom_index, pass_idx})
    return sorted(atom_set - {atom_index, pass_idx})

def bond_search(total_bond, neighbors, visited):
    found = set()
    for pair in total_bond:
        if pair[0] in neighbors and pair[1] not in visited:
            found.add(pair[1])
        elif pair[1] in neighbors and pair[0] not in visited:
            found.add(pair[0])
    if not found:
        return visited
    else:
        return bond_search(total_bond, found, visited | found)

###############################################################################


class Molecule(pybel.Molecule):
    """
    Extension of :class:`pybel.Molecule` for the generation of 3D geometries
    for structures containing more than one molecule.
    The attributes are:

    =============== ======================== ==================================
    Attribute       Type                     Description
    =============== ======================== ==================================
    `OBMol`         :class:`pybel.ob.OBMol`  An Open Babel molecule object
    `mols`          ``list``                 A list of :class:`Molecule` molecules contained in `self`
    `mols_indices`  ``list``                 Tuple of lists containing indices of atoms in the molecules
    `rotors`        ``list``                 A list of rotors
    `atom_in_rotor` ``list``                 A list to tell which rotor atoms are in
    `close_atoms`   ``list``                 A list of lists of atoms that are close together
    =============== ======================== ==================================

    Note: The molecule should have all hydrogen atoms explicitly assigned. If
    this is not the case, then segmentation faults may occur.
    """

    def __init__(self, OBMol):
        super(Molecule, self).__init__(OBMol)
        self.mols_indices = None
        self.mols = None
        self.rotors = None
        self.atom_in_rotor = None
        self.close_atoms = None

    def __getitem__(self, item):
        for atom in self:
            if item == atom.idx - 1:
                return atom
        else:
            raise IndexError('index out of range')

    def copy(self):
        """
        Create copy of `self`. The copy is somewhat reduced in that it only
        contains atoms, coords and bonds.
        """
        # Create new empty instance
        m = Molecule(pybel.ob.OBMol())
        OBMol = m.OBMol

        for atom in self:
            OBMol.AddAtom(atom.OBAtom)
        for bond in pybel.ob.OBMolBondIter(self.OBMol):
            OBMol.AddBond(bond)

        OBMol.SetTotalSpinMultiplicity(self.spin)

        m.mols_indices = self.mols_indices
        m.mols = self.mols
        return m

    def toNode(self):
        """
        Convert to :class:`node.Node` object and return the object.
        """
        atoms = []
        coords = []
        for atom in self:
            atoms.append(atom.atomicnum)
            coords.append(atom.coords)
        n = node.Node(coords, atoms, self.spin)
        n.energy = self.energy
        return n

    def toRMGSpecies(self):
        """
        Convert to :class:`rmgpy.species.Species` object and return the object.
        """
        rmg_mol = from_ob_mol(rmgpy.molecule.molecule.Molecule(), self.OBMol)
        adjlist = rmg_mol.to_adjacency_list()
        spc = Species().from_adjacency_list(adjlist)
        spc.label = ''
        return spc

    def toRMGMolecule(self):
        """
        Convert to :class:`rmgpy.molecule.Molecule` object and return the
        object.
        """
        rmg_mol = from_ob_mol(rmgpy.molecule.molecule.Molecule(), self.OBMol)
        return rmg_mol

    def assignSpinMultiplicity(self):
        """
        Assigns the spin multiplicity of all atoms based on connectivity. This
        function assumes that all hydrogens are specified explicitly.
        """
        self.OBMol.SetSpinMultiplicityAssigned()
        num_diff = 0  # Number of times that a multiplicity greater than 1 occurs
        maxspin = 1

        for atom in self:
            OBAtom = atom.OBAtom
            # Only based on difference between expected and actual valence

            #element = ELEMENT_TABLE.from_atomic_number(atom.atomicnum)
            nonbond_elec = props.valenceelec[atom.atomicnum] - \
                OBAtom.GetExplicitValence()

            # Tries to pair all electrons that are not involved in bonds, which means that singlet states are favored
            diff = nonbond_elec % 2

            if diff:
                num_diff += 1
                spin = diff + 1
                if spin > maxspin:
                    maxspin = spin
                OBAtom.SetSpinMultiplicity(spin)

        # Try to infer total spin multiplicity based on atom spins (assuming that spins in diradicals and higher are
        # opposite and favor the singlet state)
        if num_diff % 2 != 0:
            self.OBMol.SetTotalSpinMultiplicity(maxspin)
        else:
            self.OBMol.SetTotalSpinMultiplicity(1)

    def AssignSpinMultiplicity(self):
        """
        Override method in parent class.
        """
        self.assignSpinMultiplicity()

    def getH298(self, thermo_db=None):
        """
        Compute and return the standard enthalpy of formation of the structure
        in kcal/mol. A :class:`rmgpy.data.thermo.ThermoDatabase` instance can
        be supplied, which is used to search databases and use group additivity
        values.
        """
        # Load thermo database
        if thermo_db is None:
            thermo_db = ThermoDatabase()
            thermo_db.load(path.join(settings['database.directory'], 'thermo'))

        # Compute enthalpy for each molecule and add together
        H298 = 0.0
        self.separateMol()
        for mol in self.mols:
            spc = mol.toRMGSpecies()
            spc.generate_resonance_structures()
            spc.thermo = thermo_db.get_thermo_data(spc)
            H298 += spc.get_enthalpy(298.0) / _constants.KCAL2J
        # Return combined enthalpy of all molecules
        return H298

    def setCoordsFromMol(self, other):
        """
        Set the coordinates for each atom in the current molecule from the
        atoms in another one.
        """
        if len(self.atoms) != len(other.atoms):
            raise Exception('Number of atoms must match')

        for atom, other_atom in zip(self, other):
            coord_vec = other_atom.OBAtom.GetVector()
            atom.OBAtom.SetVector(coord_vec)

    def isCarbeneOrNitrene(self):
        """
        Return a boolean indicating whether or not the molecule is a carbene or
        a nitrene.
        """
        for atom in self:
            OBAtom = atom.OBAtom

            if OBAtom.GetAtomicNum() == OBElements.Carbon and OBAtom.GetExplicitValence() == 2:
                return True
            if OBAtom.GetAtomicNum() == OBElements.Nitrogen and OBAtom.GetExplicitValence() == 1:
                return True

        return False

    def optimizeGeometry(self, Qclass, **kwargs):
        """
        Perform a geometry optimization of each molecule in self using an
        electronic structure program specified in `Qclass` with the parameters
        specified in `kwargs`. Update the coordinates and energy, and return
        the number of gradient evaluations.
        """
        self.separateMol()

        if len(self.mols) > 1:
            ngrad = 0
            energy = 0.0
            name_base = kwargs.get('name', 'opt')
            err = False

            for i, mol in enumerate(self.mols):
                kwargs['name'] = name_base + str(i)
                n = mol.toNode()

                try:
                    ngrad += n.optimizeGeometry(Qclass, **kwargs)
                except QuantumError as e:
                    err, msg = True, e

                energy += n.energy
                mol_opt = n.toPybelMol()
                mol.setCoordsFromMol(mol_opt)

                self.mergeMols()

            # Raise error even if only one of the optimizations failed
            if err:
                raise QuantumError(msg)
        else:
            n = self.toNode()
            ngrad = n.optimizeGeometry(Qclass, **kwargs)
            energy = n.energy
            mol_opt = n.toPybelMol()
            self.setCoordsFromMol(mol_opt)

        self.OBMol.SetEnergy(energy)
        return ngrad

    def gen3D(self, constraint, forcefield='uff', method='ConjugateGradients', make3D=True):
        """
        Generate 3D coordinates using the specified force field.
        """
        spin = self.spin
        self.separateMol()

        # Generate 3D geometries separately
        if len(self.mols) > 1 and constraint is None:
            for mol in self.mols:
                is_hydrogen_mol = len(mol.atoms) == 2 and all(
                    a.OBAtom.GetAtomicNum() == 1 for a in mol)
                is_oxygen_mol = len(mol.atoms) == 2 and all(
                    a.OBAtom.GetAtomicNum() == 8 for a in mol)

                if make3D and len(mol.atoms) == 1:  # Atoms
                    mol.atoms[0].OBAtom.SetVector(0.0, 0.0, 0.0)
                elif make3D and is_hydrogen_mol:
                    mol.atoms[0].OBAtom.SetVector(0.0, 0.0, 0.0)
                    mol.atoms[1].OBAtom.SetVector(0.74, 0.0, 0.0)
                elif is_hydrogen_mol:
                    ind = [a.idx for a in mol]
                    bond = mol.OBMol.GetBond(ind[0], ind[1])
                    bond.SetLength(0.74)
                elif make3D and is_oxygen_mol:
                    mol.atoms[0].OBAtom.SetVector(0.0, 0.0, 0.0)
                    mol.atoms[1].OBAtom.SetVector(1.21, 0.0, 0.0)
                elif is_oxygen_mol:
                    ind = [a.idx for a in mol]
                    bond = mol.OBMol.GetBond(ind[0], ind[1])
                    bond.SetLength(1.21)
                elif constraint is None:
                    make3DandOpt(mol, forcefield=forcefield, make3D=make3D)
                else:
                    constraint_force_field(mol, constraint, forcefield, method)

            self.mergeMols()
            self.OBMol.SetTotalSpinMultiplicity(spin)
        else:
            if constraint is None:
                make3DandOpt(self, forcefield=forcefield, make3D=make3D)
            else:
                constraint_force_field(self, constraint, forcefield, method)

    def separateMol(self):
        """
        Separate molecule based on the indices in `self.mols_indices`.
        """
        if self.mols is None:
            if self.mols_indices is None:
                self.connectivityAnalysis()

            nmols = len(self.mols_indices)
            if nmols > 1:
                self.mols = []
                for mol_idx in range(nmols):
                    mol = self.copy()

                    # Obtain indices of all atoms to be deleted for current molecule and obtain corresponding atoms
                    del_indices = [atom_idx for mol_idx_2, mol_indices in enumerate(self.mols_indices)
                                   if mol_idx_2 != mol_idx
                                   for atom_idx in mol_indices]
                    del_atoms = [mol.atoms[idx].OBAtom for idx in del_indices]

                    # Delete atoms not in current molecule
                    for atom in del_atoms:
                        mol.OBMol.DeleteAtom(atom)

                    mol.assignSpinMultiplicity()  # Has to be set again because
                    # mol.OBMol.SetHydrogensAdded()
                    self.mols.append(mol)
            else:
                self.mols = [self]

    def mergeMols(self):
        """
        Merge molecules by clearing the current molecule and rewriting all
        atoms and bonds. The atoms are reordered according to the indices in
        `self.mols_indices`.
        """
        if self.mols is not None:
            spin = self.spin
            self.OBMol.Clear()

            if self.mols_indices is None:
                self.connectivityAnalysis()

            # Loop through molecules and append atoms and bonds in order
            natoms = 0
            for mol in self.mols:
                for atom in mol:
                    self.OBMol.AddAtom(atom.OBAtom)
                for bond in pybel.ob.OBMolBondIter(mol.OBMol):
                    self.OBMol.AddBond(
                        bond.GetBeginAtomIdx() + natoms, bond.GetEndAtomIdx() +
                        natoms, bond.GetBondOrder()
                    )
                natoms += len(mol.atoms)

            # Reorder atoms
            mols_indices_new = [atom_idx for mol_indices in self.mols_indices for atom_idx in mol_indices]
            neworder = natoms * [0]
            for i, atom_idx in enumerate(mols_indices_new):
                neworder[atom_idx] = i + 1
            self.OBMol.RenumberAtoms(neworder)

            # self.OBMol.SetHydrogensAdded()
            self.OBMol.SetTotalSpinMultiplicity(spin)

    def connectivityAnalysis(self):
        """
        Analyze bonds to determine which atoms are connected and form a
        molecule.
        """
        # Extract bonds
        bonds = [[bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1]
                 for bond in pybel.ob.OBMolBondIter(self.OBMol)]
        if bonds:
            # Create first molecular fragment from first bond and start keeping track of atoms
            molecules = [bonds[0][:]]
            atoms_used = bonds[0][:]

            # Loop over remaining bonds
            for bond in bonds[1:]:
                ind1, ind2 = -1, -2

                for idx, molecule in enumerate(molecules):
                    if bond[0] in molecule:
                        ind1 = idx
                    if bond[1] in molecule:
                        ind2 = idx

                # Skip bond if both atoms are already contained in the same molecule
                if ind1 == ind2:
                    continue
                # Combine fragments if they are connected through bond
                if ind1 != -1 and ind2 != -2:
                    molecules[ind1].extend(molecules[ind2])
                    del molecules[ind2]
                # Add new atom to fragment if it is connected
                elif ind1 != -1:
                    molecules[ind1].append(bond[1])
                    atoms_used.append(bond[1])
                elif ind2 != -2:
                    molecules[ind2].append(bond[0])
                    atoms_used.append(bond[0])
                # Add new fragment if it does not connect to any other ones
                else:
                    molecules.append(bond)
                    atoms_used.extend(bond)

            # Add atoms that are not involved in bonds
            for atom in range(len(self.atoms)):
                if atom not in atoms_used:
                    molecules.append([atom])

            # Sort molecules and store result
            self.mols_indices = tuple(sorted(molecule) for molecule in molecules)
        else:
            self.mols_indices = tuple([atom] for atom in range(len(self.atoms)))

    def detRotors(self, mols_indices, fixed_atoms=[]):
        """
        Determine the rotors and atoms in the rotors of the molecule.
        """
        self.rotors = []
        self.atom_in_rotor = []
        natoms = len(self.atoms)
        total_bond = [(bond_1.GetBeginAtomIdx() - 1, bond_1.GetEndAtomIdx() - 1) for bond_1 in pybel.ob.OBMolBondIter(self.OBMol)]

        for bond_1 in pybel.ob.OBMolBondIter(self.OBMol):
            if bond_1.IsRotor():
                ref_1, ref_2 = bond_1.GetBeginAtomIdx() - 1, bond_1.GetEndAtomIdx() - 1
                if mols_indices[ref_1] not in fixed_atoms or mols_indices[ref_2] not in fixed_atoms:
                    self.rotors.append((ref_1, ref_2))
                    atom_in_rotor = [False] * natoms

                    conns = det_connectivity(total_bond, atom_index=ref_1, pass_idx=ref_2)
                    for conn in conns:
                        if conn in fixed_atoms:
                            atom_in_rotor[ref_2] = True
                            break
                    else:
                        atom_in_rotor[ref_1] = True

                    new_atom = True
                    while new_atom:
                        new_atom = False
                        for bond_2 in pybel.ob.OBMolBondIter(self.OBMol):
                            ref_3, ref_4 = bond_2.GetBeginAtomIdx() - 1, bond_2.GetEndAtomIdx() - 1
                            if ref_1 == ref_3 and ref_2 == ref_4:
                                continue
                            if ref_2 == ref_3 and ref_1 == ref_4:
                                continue
                            if atom_in_rotor[ref_3] ^ atom_in_rotor[ref_4]:
                                atom_in_rotor[ref_3], atom_in_rotor[ref_4] = True, True
                                new_atom = True

                    self.atom_in_rotor.append(atom_in_rotor)


    def detCloseAtoms(self, d):
        """
        Find atoms that are closer together than `d`.
        """
        natoms = len(self.atoms)
        self.close_atoms = [[0 for i in range(natoms)] for k in range(natoms)]
        for i in range(0, natoms-1):
            for j in range(i + 1, natoms):
                if self.atoms[i].OBAtom.GetDistance(self.atoms[j].OBAtom) < d:
                    self.close_atoms[i][j], self.close_atoms[j][i] = 1, 1

###############################################################################


class Arrange3D(object):
    """
    Arranging of :class:`Molecule` or :class:`pybel.Molecule` molecule objects in
    3D space by global translations and rotations and internal rotations.
    The attributes are:

    ============== ================== =========================================
    Attribute      Type               Description
    ============== ================== =========================================
    `mol_1`        :class:`Molecule`  The first molecule
    `mol_2`        :class:`Molecule`  The second molecule
    `bonds_1`      ``list``           The bonds in molecule 1
    `bonds_2`      ``list``           The bonds in molecule 2
    `torsions_1`   ``list``           The torsions in molecule 1
    `torsions_2`   ``list``           The torsions in molecule 2
    `d_intermol`   ``float``          The intermolecular distance
    `d_intramol`   ``float``          The intramolecular distance
    `nodes_1`      ``list``           List of :class:`node.Node` in molecule 1
    `nodes_2`      ``list``           List of :class:`node.Node` in molecule 2
    ============== ================== =========================================

    """

    def __init__(self, mol_1, mol_2, fixed_atoms=None, cluster_bond = None):
        # set initial position has been disable so i think we can remove this fragment constraint
        # if not (0 < len(mol_1.mols) <= 4 and 0 < len(mol_2.mols) <= 4):
        #raise Exception('More than 4 molecules are not supported')

        self.mol_1 = None
        self.mol_2 = None
        self.bonds_1 = None
        self.bonds_2 = None
        self.torsions_1 = None
        self.torsions_2 = None
        self.d_intermol = None
        self.d_intramol = None
        self.dof_1 = None
        self.def_2 = None
        self.nodes_1 = None
        self.nodes_2 = None
        self.fdist_1 = None
        self.fdist_2 = None

        if fixed_atoms is None:
            self.fixed_atoms = []
        else:
            self.fixed_atoms = fixed_atoms

        if cluster_bond is None:
            self.cluster_bond = []
        else:
            self.cluster_bond = cluster_bond
        self.initializeVars(mol_1, mol_2)

    def initializeVars(self, mol_1, mol_2, d_intermol=3.0, d_intramol=2.0):
        """
        Set up class variables and determine the bonds and torsions to be
        matched between reactant and product.
        """
        self.bonds_1, self.bonds_2, self.torsions_1, self.torsions_2 = [], [], [], []
        self.mol_1, self.mol_2 = mol_1, mol_2
        self.d_intermol = d_intermol
        self.d_intramol = d_intramol

        bonds_mol_1 = [tuple(sorted((bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1)))
                       for bond in pybel.ob.OBMolBondIter(mol_1.OBMol)]

        bonds_mol_2 = [tuple(sorted((bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1)))
                       for bond in pybel.ob.OBMolBondIter(mol_2.OBMol)]

        # atom_in_mol tells which mol the atom is in
        atom_in_mol_1 = [0] * len(mol_1.atoms)
        atom_in_mol_2 = [0] * len(mol_2.atoms)

        for i, mol_indices in enumerate(mol_1.mols_indices):
            for j in mol_indices:
                atom_in_mol_1[j] = i

        for i, mol_indices in enumerate(mol_2.mols_indices):
            for j in mol_indices:
                atom_in_mol_2[j] = i

        # Initialize broken_bonds list
        for bond in list(set(bonds_mol_2) ^ set(bonds_mol_1)):
            i, j = atom_in_mol_1[bond[0]], atom_in_mol_1[bond[1]]
            self.bonds_1.append([(i, mol_1.mols_indices[i].index(bond[0])), (j, mol_1.mols_indices[j].index(bond[1]))])
            i, j = atom_in_mol_2[bond[0]], atom_in_mol_2[bond[1]]
            self.bonds_2.append([(i, mol_2.mols_indices[i].index(bond[0])), (j, mol_2.mols_indices[j].index(bond[1]))])

        # Connectivity tells whether two atoms are linked by a bond
        natoms = len(mol_1.atoms)
        connectivity = [[0 for i in range(natoms)] for k in range(natoms)]

        for bond in list(set(bonds_mol_1) | set(bonds_mol_2)):
            connectivity[bond[0]][bond[1]], connectivity[bond[1]][bond[0]] = 1, 1

        # Generate a list of dihedral angles to be compared between reactant and product
        for i in range(natoms):
            for j in range(natoms):
                if not connectivity[i][j]:
                    continue
                for k in range(natoms):
                    if (not connectivity[j][k]) or k == i:
                        continue
                    for l in range(natoms):
                        if (not connectivity[k][l]) or l == j or l <= i:
                            continue
                        a, b, c, d = atom_in_mol_1[i], atom_in_mol_1[j], atom_in_mol_1[k], atom_in_mol_1[l]
                        # if i not in self.constraint and j not in self.constraint and k not in self.constraint and l not in self.constraint:  # debug
                        self.torsions_1.append([
                            (a, mol_1.mols_indices[a].index(i)),
                            (b, mol_1.mols_indices[b].index(j)),
                            (c, mol_1.mols_indices[c].index(k)),
                            (d, mol_1.mols_indices[d].index(l))
                        ])
                        a, b, c, d = atom_in_mol_2[i], atom_in_mol_2[j], atom_in_mol_2[k], atom_in_mol_2[l]
                        self.torsions_2.append([
                            (a, mol_2.mols_indices[a].index(i)),
                            (b, mol_2.mols_indices[b].index(j)),
                            (c, mol_2.mols_indices[c].index(k)),
                            (d, mol_2.mols_indices[d].index(l))
                        ])

        # Determine rotors and connectivity in molecules
        self.dof_1, self.def_2 = 6 * (len(self.mol_1.mols) - 1), 6 * (len(self.mol_2.mols) - 1)

        for idx, mol in enumerate(mol_1.mols):
            mol.detRotors(mol.mols_indices[idx], self.fixed_atoms)
            mol.detCloseAtoms(d_intramol)
            self.dof_1 += len(mol.rotors)

        for idx, mol in enumerate(mol_2.mols):
            mol.detRotors(mol.mols_indices[idx], self.fixed_atoms)
            mol.detCloseAtoms(d_intramol)
            self.def_2 += len(mol.rotors)

        # Convert mols to nodes and center molecules
        self.nodes_1 = [mol.toNode() for mol in self.mol_1.mols]
        self.nodes_2 = [mol.toNode() for mol in self.mol_2.mols]
        self.nodes_1_copy = self.nodes_1[:]
        self.nodes_2_copy = self.nodes_2[:]
        if not self.fixed_atoms:
            self.setInitialPositions(self.nodes_1)
            self.setInitialPositions(self.nodes_2)

        if len(self.nodes_1) > 1:
            fd1 = [node.getCentroid() for node in self.nodes_1]

            #center = [node.getCentroid() for node in self.nodes_1]
            #center = [node.getCenterOfMass() for node in self.nodes_1]
            self.fdist_1 = [np.linalg.norm(a-b) for a, b in zip(fd1, fd1[1:] + fd1[:-1])]
        elif len(self.nodes_2) > 1:
            fd2 = [node.getCentroid() for node in self.nodes_1]
            self.fdist_2 = [np.linalg.norm(a-b) for a, b in zip(fd2, fd2[1:] + fd2[:-1])]

    def arrangeIn3D(self):
        """
        Arrange the molecules in 3D-space by global translations and rotations
        and internal rotor rotations. Return an empty string if optimization
        converged or if there are no degrees of freedom for optimization, the
        status message otherwise
        """
        ret = ''
        dof = self.dof_1 + self.def_2
        if dof != 0:
            #a = self.mol_1
            #b = self.nodes_1
            """
            def callbackF(Xi):
                print(self.objectiveFunction(Xi[:dof]))
                coords_1 = self.newCoords(self.mol_1.mols, self.nodes_1, Xi[:self.dof_1])
                with open('test.xyz', 'a') as f:
                    for i in range(0, 2):
                        b[i].coords = coords_1[i]
                        a.mols[i].setCoordsFromMol(b[i].toPybelMol())
                        f.write('{}\n\n{}\n'.format(str(36), str(a.toNode())))
            """

            disps_guess = np.array([0.0]*dof)
            # The 10 x 10 box should find all minimum
            #bounds = Bounds([-5.0]*dof, [5.0]*dof)
            if self.fixed_atoms == []:
                result = minimize(self.objectiveFunction, disps_guess,
                                        constraints={'type': 'ineq', 'fun': self.constraintFunction},
                                        method='SLSQP',
                                        options={'disp': False, 'maxiter':1000})  # , callback = callbackF, 'eps':1e-10
            else:
                # result = minimize(self.objectiveFunction, disps_guess,
                #                         constraints=[{'type': 'ineq', 'fun': self.constraintFunction},
                #                                     {'type': 'ineq', 'fun': self.third_constraintFunction}
                #                                     ],
                #                         method='COBYLA',
                #                         options={'disp': False})

                result = minimize(self.objectiveFunction, disps_guess,
                                        constraints=[{'type': 'ineq', 'fun': self.constraintFunction},
                                                    {'type': 'ineq', 'fun': self.third_constraintFunction}
                                                    ],
                                        method='SLSQP',
                                        options={'disp': False, 'maxiter':200, 'ftol':0.01})  # , callback = callbackF, 'eps':1e-10

                # if not s_result.success:
                #     result = result
                # else:
                #     result = s_result

            assert result.success

            if not result.success:
                message = ('Optimization in arrangeIn3D terminated with status ' +
                           str(result.status) + ':\n' + result.message + '\n')
                ret = message

            coords_1 = self.newCoords(self.mol_1.mols, self.nodes_1, result.x[:self.dof_1])
            coords_2 = self.newCoords(self.mol_2.mols, self.nodes_2, result.x[self.dof_1:])

            # coords_1 = self.backtransform(coords_1, self.nodes_1_copy)
            # coords_2 = self.backtransform(coords_2, self.nodes_2_copy)

            for i in range(0, len(self.mol_1.mols)):
                self.nodes_1[i].coords = coords_1[i]
                self.mol_1.mols[i].setCoordsFromMol(self.nodes_1[i].toPybelMol())

            for i in range(0, len(self.mol_2.mols)):
                self.nodes_2[i].coords = coords_2[i]
                self.mol_2.mols[i].setCoordsFromMol(self.nodes_2[i].toPybelMol())

            if len(self.mol_1.mols) > 1:
                self.mol_1.mergeMols()
            if len(self.mol_2.mols) > 1:
                self.mol_2.mergeMols()

        return ret

    @staticmethod
    def centerAndFindDistances(nodes):
        """
        Center the molecules about the origin and return the distances between
        the origin and the atom farthest from the origin, which can be used as
        size estimates for the molecules.
        """
        max_distances = []
        for n in nodes:
            disp = n.getCentroid()
            n.translate(-disp)
            coords = n.coords
            max_distance = 0.0
            for coord in coords:
                distance = np.linalg.norm(coord)
                if distance > max_distance:
                    max_distance = distance
            max_distances.append(max_distance)
        return max_distances

    @staticmethod
    def fragment_dist(nodes):
        """
        Center the molecules about the origin and return the distances between
        the origin and the atom farthest from the origin, which can be used as
        size estimates for the molecules.
        """
        disp = nodes.getCentroid()
        return disp

    def setInitialPositions(self, nodes):
        """
        Determine initial positions of molecules
        """
        sizes = self.centerAndFindDistances(nodes)
        d = self.d_intermol
        nmols = len(nodes)
        disp = []
        if nmols == 2:
            t = np.array([d + sizes[0] + sizes[1], 0.0, 0.0])
            disp.append(t)
        # Arrange molecules in triangle if there are three molecules
        elif nmols == 3:
            t1 = np.array([d + sizes[0] + sizes[1], 0.0, 0.0])
            d1 = d + sizes[0] + sizes[1]
            d2 = d + sizes[0] + sizes[2]
            d3 = d + sizes[1] + sizes[2]
            y = (-d1 ** 4.0 + 2.0 * d1 ** 2.0 * d2 ** 2.0 + 2.0 * d1 ** 2.0 * d3 ** 2.0 -
                 d2 ** 4.0 + 2.0 * d2 ** 2.0 * d3 ** 2.0 - d3 ** 4.0) ** 0.5 / (2.0 * d1)
            x = (d2 ** 2.0 - y ** 2.0) ** 0.5
            t2 = np.array([x, -y, 0.0])
            disp += [t1] + [t2]
        # Arrange molecules in square if there are four molecules
        elif nmols == 4:
            x = max(d + sizes[0] + sizes[1], d + sizes[2] + sizes[3])
            y = max(d + sizes[0] + sizes[2], d + sizes[1] + sizes[3])
            t1 = np.array([x, 0.0, 0.0])
            t2 = np.array([0.0, -y, 0.0])
            t3 = np.array([x, -y, 0.0])
            disp += [t1] + [t2] + [t3]

        for i in range(1, nmols):
            nodes[i].translate(disp[i - 1])

    @staticmethod
    def calcBondLens(coords, bonds):
        """
        Calculate the lengths of the bonds in a given structure.
        """
        dist = []
        for bond in bonds:
            coord_vect_1 = coords[bond[0][0]][bond[0][1]]
            coord_vect_2 = coords[bond[1][0]][bond[1][1]]
            diff = coord_vect_1 - coord_vect_2
            dist.append(np.linalg.norm(diff))
        return dist

    @staticmethod
    def calcDihedralAngs(coords, torsions):
        """
        Calculate the dihedral angles of the torsions in a given structure.
        """
        angles = []
        for d in torsions:
            a1 = coords[d[0][0]][d[0][1]]
            a2 = coords[d[1][0]][d[1][1]]
            a3 = coords[d[2][0]][d[2][1]]
            a4 = coords[d[3][0]][d[3][1]]
            b1 = a1 - a2
            b2 = a2 - a3
            b3 = a3 - a4
            n1 = np.cross(b1, b2)
            n1 /= np.linalg.norm(n1)
            n2 = np.cross(b2, b3)
            n2 /= np.linalg.norm(n2)
            m1 = np.cross(n1, b2)
            m1 /= np.linalg.norm(m1)
            x = np.dot(n1, n2)
            y = np.dot(m1, n2)
            angles.append(np.arctan2(x, y))
        return angles

    @staticmethod
    def minIntermolDist(coords):
        """
        Determine minimum distance between molecules.
        """
        nmols = len(coords)
        dist_min = 0

        for i in range(0, nmols - 1):
            for j in range(i + 1, nmols):
                for k, coord_vect_1 in enumerate(coords[i]):
                    for l, coord_vect_2 in enumerate(coords[j]):
                        diff = coord_vect_1 - coord_vect_2
                        dist = np.linalg.norm(diff)
                        if dist_min > dist or dist_min == 0:
                            dist_min = dist
        return dist_min

    @staticmethod
    def minIntramolDist(coords, mols):
        """
        Determine minimum distance between two nonbonded atoms in a molecule.
        """
        dist_min = 0.0
        for i, mol in enumerate(mols):
            natoms = len(coords[i])
            for j in range(0, natoms - 1):
                for k in range(j + 1, natoms):
                    if not mol.close_atoms[j][k]:
                        diff = coords[i][j] - coords[i][k]
                        dist = np.linalg.norm(diff)
                        if dist_min > dist or dist_min == 0:
                            dist_min = dist
        return dist_min

    @staticmethod
    def translate(coords, trans_vec):
        """
        Translate all atoms in the molecular configuration by `trans_vec`,
        which is of type :class:`numpy.ndarray` and of size 3 x 1.
        """
        coords = coords + trans_vec
        return coords

    @staticmethod
    def rotate(coords, rot_mat):
        """
        Rotate molecular configuration about the origin using orthogonal
        rotation matrix `rot_mat` which is of type :class:`numpy.ndarray`
        and of size 3 x 3.
        """
        coords = (rot_mat.dot(coords.T)).T
        return coords

    @staticmethod
    def get_idx(target, array):
        for i, mol_fragment in enumerate(array):
            try:
                j = mol_fragment.index(target)
            except ValueError:
                continue
            yield i, j

    def rotateRotor(self, coords, angle, rotor, atom_in_rotor):
        """
        Rotate internal rotor. Note that global rotation is inevitable.
        """
        n = len(coords)
        centroid = coords.sum(axis=0) / n
        disp = coords[rotor[0]]
        coords = self.translate(coords, -disp)
        axis = coords[rotor[0]] - coords[rotor[1]]
        rot_mat = util.rotationMatrix(angle, axis)
        for j in range(n):
            if atom_in_rotor[j]:
                coords[j] = self.rotate(coords[j], rot_mat)
        coords = self.translate(coords, +disp)
        return coords

    def rotateMol(self, coords, angles):
        """
        Rotate molecular configuration.
        """
        n, m = coords.shape
        centroid = coords.sum(axis=0) / n
        rot_mat = util.rotationMatrix(angles)
        coords = self.translate(coords, -centroid)
        coords = self.rotate(coords, rot_mat)
        coords = self.translate(coords, centroid)
        return coords

    def newCoords(self, mols, nodes, disps):
        """
        Determine new coordinates after internal rotor rotations and global
        rotations and translations.
        """
        nmols = len(mols)
        trans_disps = disps[:3 * (nmols - 1)]
        rot_disps = disps[3 * (nmols - 1):6 * (nmols - 1)]
        tort_disps = disps[6 * (nmols - 1):]
        coords_all = []
        nrots = 0
        
        if nmols == 1:
            mol = mols[0]
            coords = nodes[0].coords
            for j in range(len(mol.rotors)):
                coords = self.rotateRotor(coords, tort_disps[nrots], mol.rotors[j], mol.atom_in_rotor[j])
                nrots += 1
            coords_all.append(coords)
        else:
            check = []
            for i in range(0, nmols):
                index = mols[i].mols_indices[i]
                for idx in index:
                    if idx in self.fixed_atoms:
                        check.append(i)
            check = set(check)

            for i in range(0, nmols):
                mol = mols[i]
                coords = nodes[i].coords
                for j in range(len(mol.rotors)):
                    coords = self.rotateRotor(coords, tort_disps[nrots], mol.rotors[j], mol.atom_in_rotor[j])
                    nrots += 1

                if i not in check:
                    if i > 0:
                        coords = self.rotateMol(coords, rot_disps[3 * (i - 1):3 * i])
                        coords = self.translate(coords, trans_disps[3 * (i - 1):3 * i])
                    else:
                        coords = self.rotateMol(coords, rot_disps[3 * i:3 * (i + 1)])
                        coords = self.translate(coords, trans_disps[3 * i:3 * (i + 1)])              
                coords_all.append(coords)

        return coords_all

    def gradient(self, disps):
        """
        Gradient for scipy minimize. (jac)
        """
        coords_1 = self.newCoords(self.mol_1.mols, self.nodes_1, disps[:self.dof_1])
        coords_2 = self.newCoords(self.mol_2.mols, self.nodes_2, disps[self.dof_1:])
        atom_symbol_1, atom_symbol_2 = [], []
        for mol in self.mol_1.mols:
            atoms = tuple(atom.atomicnum for atom in mol)
            for atom in atoms:
                atom_symbol_1.append(ob.GetSymbol(atom))
        for mol in self.mol_2.mols:
            atoms = tuple(atom.atomicnum for atom in mol)
            for atom in atoms:
                atom_symbol_2.append(ob.GetSymbol(atom))
        force_1 = get_ob_forces(atom_symbol_1, np.concatenate(coords_1), cluster_bond = self.cluster_bond, freeze_index=self.fixed_atoms)
        force_2 = get_ob_forces(atom_symbol_2, np.concatenate(coords_2), cluster_bond = self.cluster_bond, freeze_index=self.fixed_atoms)

    def backtransform(self, new_coords, nodes_copy):
        """
        Let the constrained atom go back to intial prosition after arrange
        """
        matches = []
        for idx, i in enumerate(self.mol_1.mols_indices):
            for constrained_atom in self.fixed_atoms:
                if constrained_atom in i:
                    target1 = i.index(constrained_atom)
                    matches.append((idx, target1))
        coords_1, coords_2 = [], []

        for match in matches:
            a = nodes_copy[match[0]].coords[match[1]]
            b = new_coords[match[0]][match[1]]
            coords_1.append(a)
            coords_2.append(b)

        center_1 = np.array(coords_1).sum(axis=0) / len(matches) # initial
        center_2 = np.array(coords_2).sum(axis=0) / len(matches) # after arrange
        disp = center_2 - center_1
        coords_all = []
        for coord in new_coords:
            coords = self.translate(coord, -disp)
            coords_all.append(coords)

        return coords_all

    def objectiveFunction(self, disps):
        """
        Sum of the difference of the bond lengths and dihedral angles between
        reactant and product.
        """
        coords_1 = self.newCoords(self.mol_1.mols, self.nodes_1, disps[:self.dof_1])
        coords_2 = self.newCoords(self.mol_2.mols, self.nodes_2, disps[self.dof_1:])
        b1 = self.calcBondLens(coords_1, self.bonds_1)
        b2 = self.calcBondLens(coords_2, self.bonds_2)
        d1 = self.calcDihedralAngs(coords_1, self.torsions_1)
        d2 = self.calcDihedralAngs(coords_2, self.torsions_2)

        # atom_symbol_1, atom_symbol_2 = [], []
        # for mol in self.mol_1.mols:
        #     atoms = tuple(atom.atomicnum for atom in mol)
        #     for atom in atoms:
        #         atom_symbol_1.append(ob.GetSymbol(atom))
        # for mol in self.mol_2.mols:
        #     atoms = tuple(atom.atomicnum for atom in mol)
        #     for atom in atoms:
        #         atom_symbol_2.append(ob.GetSymbol(atom))
        # force_1 = get_ob_forces(atom_symbol_1, np.concatenate(coords_1), cluster_bond = self.cluster_bond, freeze_index=self.fixed_atoms)
        # force_2 = get_ob_forces(atom_symbol_2, np.concatenate(coords_2), cluster_bond = self.cluster_bond, freeze_index=self.fixed_atoms)

        # with open ('./visual_arrange.xyz', 'a') as f:
        #     coords = "\n".join(["{} {:10.08f} {:10.08f} {:10.08f}".format(a, *c) for a, c in zip(atom_symbol_1, np.concatenate(coords_1))])
        #     f.write(f"{len(atom_symbol_1)}\n\n{coords}\n")

        val_b, val_d = 0.0, 0.0
        for i in range(len(b1)):
            val_b += np.abs(b1[i]-b2[i])
        for i in range(len(d1)):
            d = d1[i] - d2[i]
            if d > np.pi:
                val_d += np.abs(d - 2 * np.pi)
            elif d < -np.pi:
                val_d += np.abs(d + 2 * np.pi)
            else:
                val_d += np.abs(d)

        val = 5 * val_d + val_b #+ (force_1 + force_2) * 1e-6

        return val

    def constraintFunction(self, disps):
        """
        Optimization constraint function to make sure the shortest
        intermolecular and intramolecular distances between two atoms are
        greater than user specified values.
        """
        coords_1 = self.newCoords(self.mol_1.mols, self.nodes_1, disps[:self.dof_1])
        coords_2 = self.newCoords(self.mol_2.mols, self.nodes_2, disps[self.dof_1:])

        intermol_dists = [self.minIntermolDist(coords_1), self.minIntermolDist(coords_2)]
        intramol_dists = [self.minIntramolDist(coords_1, self.mol_1.mols),
                          self.minIntramolDist(coords_2, self.mol_2.mols)]

        val = min([a - self.d_intermol for a in intermol_dists if a != 0] +
                  [b - self.d_intramol for b in intramol_dists if b != 0])

        return val *.1

    def second_constraintFunction(self, disps):
        """
        The distance between constrained atoms should be the same after arranging
        """
        dist1, dist2 = 0.0, 0.0

        coords_1 = self.newCoords(self.mol_1.mols, self.nodes_1, disps[:self.dof_1])
        coords_2 = self.newCoords(self.mol_2.mols, self.nodes_2, disps[self.dof_1:])

        # combs = combinations(self.fixed_atoms, 2)
        for constrained_atom in self.fixed_atoms:
            matches_1 = [match for match in self.get_idx(constrained_atom, self.mol_1.mols_indices)]
            matches_2 = [match for match in self.get_idx(constrained_atom, self.mol_2.mols_indices)]
            dist1 += np.linalg.norm(coords_1[matches_1[0][0]][matches_1[0][1]] - self.nodes_1_copy[matches_1[0][0]].coords[matches_1[0][1]])
            dist2 += np.linalg.norm(coords_2[matches_2[0][0]][matches_2[0][1]] - self.nodes_2_copy[matches_2[0][0]].coords[matches_2[0][1]])
        return (dist1 + dist2) * 1e-3

    def third_constraintFunction(self, disps, threshold = 0.0):
        """
        Make sure the fragments without bonding between other will not go far away during arranging.
        """
        val_dist1, val_dist2 = 0.0, 0.0
        coords_1 = self.newCoords(
            self.mol_1.mols, self.nodes_1, disps[:self.dof_1])
        coords_2 = self.newCoords(
            self.mol_2.mols, self.nodes_2, disps[self.dof_1:])
        fragment_dist_1 = [node.sum(axis=0)/len(node) for node in coords_1]
        fragment_dist_2 = [node.sum(axis=0)/len(node) for node in coords_2]
        a = [np.linalg.norm(a-b) for a, b in zip(fragment_dist_1,
                                                 fragment_dist_1[1:] + fragment_dist_1[:-1])]
        b = [np.linalg.norm(a-b) for a, b in zip(fragment_dist_2,
                                                 fragment_dist_2[1:] + fragment_dist_2[:-1])]

        # if len(self.nodes_1) > 1 and len(self.nodes_2) > 1:
        #     threshold *= 2

        # if len(self.nodes_1) > 1:
        #     for i in range(len(a)):
        #         val_dist1 += abs(self.fdist_1[i] - a[i])
        #     val_dist1 /=len(a)
        # elif len(self.nodes_2) > 1:
        #     for i in range(len(b)):
        #         val_dist2 += abs(self.fdist_2[i] - b[i])
        #     val_dist2 /= len(b)

        # return val_dist1 + val_dist2

        if len(self.nodes_1) > 1:
            for i in range(len(a)):
                val_dist1 += abs(self.fdist_1[i] - a[i])
            val_dist1 /=len(a)
        elif len(self.nodes_2) > 1:
            for i in range(len(b)):
                val_dist2 += abs(self.fdist_2[i] - b[i])
            val_dist2 /= len(b)
        return - (val_dist1 + val_dist2 - threshold)

    def second_constraintFunction_cobyla_1(self, disps):
        """
        The distance between constrained atoms should be the same after arranging
        """
        dist1, dist2 = 0.0, 0.0

        coords_1 = self.newCoords(self.mol_1.mols, self.nodes_1, disps[:self.dof_1])
        coords_2 = self.newCoords(self.mol_2.mols, self.nodes_2, disps[self.dof_1:])

        # combs = combinations(self.fixed_atoms, 2)
        for constrained_atom in self.fixed_atoms:
            matches_1 = [match for match in self.get_idx(constrained_atom, self.mol_1.mols_indices)]
            matches_2 = [match for match in self.get_idx(constrained_atom, self.mol_2.mols_indices)]
            dist1 += np.linalg.norm(coords_1[matches_1[0][0]][matches_1[0][1]] - self.nodes_1_copy[matches_1[0][0]].coords[matches_1[0][1]])
            dist2 += np.linalg.norm(coords_2[matches_2[0][0]][matches_2[0][1]] - self.nodes_2_copy[matches_2[0][0]].coords[matches_2[0][1]])
        return dist1 + dist2 - 1e-6

    def second_constraintFunction_cobyla_2(self, disps):
        """
        The distance between constrained atoms should be the same after arranging
        """
        dist1, dist2 = 0.0, 0.0

        coords_1 = self.newCoords(self.mol_1.mols, self.nodes_1, disps[:self.dof_1])
        coords_2 = self.newCoords(self.mol_2.mols, self.nodes_2, disps[self.dof_1:])

        # combs = combinations(self.fixed_atoms, 2)
        for constrained_atom in self.fixed_atoms:
            matches_1 = [match for match in self.get_idx(constrained_atom, self.mol_1.mols_indices)]
            matches_2 = [match for match in self.get_idx(constrained_atom, self.mol_2.mols_indices)]
            dist1 += np.linalg.norm(coords_1[matches_1[0][0]][matches_1[0][1]] - self.nodes_1_copy[matches_1[0][0]].coords[matches_1[0][1]])
            dist2 += np.linalg.norm(coords_2[matches_2[0][0]][matches_2[0][1]] - self.nodes_2_copy[matches_2[0][0]].coords[matches_2[0][1]])
        return - (dist1 + dist2) + 1e-6