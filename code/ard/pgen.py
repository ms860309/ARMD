#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contains the :class:`Generate` for generating product structures.
"""

# buildin
from operator import itemgetter
from itertools import combinations, count

# third party
from openbabel import pybel
import networkx as nx

# local application imports
import props
import _constants
import gen3D
import numpy as np
from graph import make_graph, is_isomorphic
from bond_rearrangement import get_bond_rearrangs, BondRearrangement
from atoms import Atom
from species import Species
from dataclasses import dataclass

ELEMENT_TABLE = props.ElementData()
###############################################################################


class StructureError(Exception):
    """
    An exception class for errors that occur while generating structures.
    """
    pass

###############################################################################


class Generate(object):
    """
    Generation of product structures.
    The attributes are:

    ===============         ========================== ================================
    Attribute               Type                       Description
    ===============         ========================== ================================
    `reac_mol`              class:`gen3D.Molecule`     A molecule object for the reactant structure
    `reactant_inchikey`     ``str``                    A string of reactant inchi key is use to prevent product equal to reactant
    `reac_mol_graph`        ``mol_graph.graph``        A graph contain bond connectivity and atoms are use to filter the invalid products
    `bond_dissociation_cutoff` ``float``               Use to approximate barrier (Over estimate).
    `atoms`                 ``tuple``                  A tuple containing the atomic numbers of reactant/product structures
    `prod_mols`             ``list``                   A list of :class:`gen3D.Molecule` product structures
    `add_bonds`             ``list(tuple)``            Is for SSM
    `break_bonds`           ``list(tuple)``            Is for SSM
    `use_inchi_key`         ``int``                    Prevent the products equal to reactant, the graph is not always filter them
    
    =============== ========================== ================================

    Note: Bonds are represented as
          (beginAtomIdx, endAtomIdx, bondOrder)
    """

    def __init__(self, reac_mol, fixed_atom=None, **kwargs):
        self.reac_mol = reac_mol
        self.reactant_inchikey = [reac_mol.write('inchiKey').strip()]
        self.reac_mol_graph = kwargs['graph']
        self.bond_dissociation_cutoff = float(kwargs['bond_dissociation_cutoff'])
        self.nform = kwargs['nform']
        self.nbreak = kwargs['nbreak']
        self.atoms = None
        self.prod_mols, self.add_bonds, self.break_bonds = [], [], []
        if kwargs['use_inchi_key'] == '1':
            self.use_inchi_key = True
        else:
            self.use_inchi_key = False
        self.fixed_atom = fixed_atom

        self.initialize()
        self.reactant_coords = [np.array([atom.coords for atom in self.reac_mol]).reshape(len(self.atoms), 3)]
        self.catalyst = kwargs['catalyst']

    def initialize(self):
        """
        Set the canonical SMILES for the reactant and extract the atomic
        numbers.
        """
        self.atoms = tuple(atom.atomicnum for atom in self.reac_mol)
        reactant_bonds = [tuple(sorted((bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1)) + [bond.GetBondOrder()])
                            for bond in pybel.ob.OBMolBondIter(self.reac_mol.OBMol)]
        self.reactant_bonds = tuple(sorted(reactant_bonds))

        if self.fixed_atom:
            self.active_site_oxygen = [active_site_atom for active_site_atom in self.fixed_atom if self.atoms[active_site_atom] == 8]
            self.active_site_metal = [active_site_atom for active_site_atom in self.fixed_atom if self.atoms[active_site_atom] in [42, 50, 74, 13]]
            self.reactant_oxygen = [idx for idx, reactant_atoms in enumerate(self.atoms) if idx not in self.fixed_atom and reactant_atoms == 8]
            self.reactant_carbon = [idx for idx, reactant_atoms in enumerate(self.atoms) if idx not in self.fixed_atom and reactant_atoms == 6]
            # The zeolite usually reactant with its faced active site oxygen
            self.active_site_oxygen.remove(12)
            # Extract bonds as an unmutable sequence (indices are made compatible with atom list)
            # self.reactant_hydrogen is -OH self.reactant_hydrogen_2 is C-H
            reactant_hydrogen_1 =[bond[0] for bond in self.reactant_bonds 
                                if bond[0] not in self.fixed_atom and bond[1] not in self.fixed_atom and self.atoms[bond[0]] == 1 and self.atoms[bond[1]] == 8]
            self.reactant_hydrogen =[bond[1] for bond in self.reactant_bonds 
                                if bond[0] not in self.fixed_atom and bond[1] not in self.fixed_atom and self.atoms[bond[1]] == 1 and self.atoms[bond[0]] == 8]
            self.reactant_hydrogen.extend(reactant_hydrogen_1)
            reactant_hydrogen_1 =[bond[0] for bond in self.reactant_bonds 
                                if bond[0] not in self.fixed_atom and bond[1] not in self.fixed_atom and self.atoms[bond[0]] == 1 and self.atoms[bond[1]] == 6]
            self.reactant_hydrogen_2 =[bond[1] for bond in self.reactant_bonds 
                                if bond[0] not in self.fixed_atom and bond[1] not in self.fixed_atom and self.atoms[bond[1]] == 1 and self.atoms[bond[0]] == 6]
            self.reactant_hydrogen_2.extend(reactant_hydrogen_1)
            active_site_hydrogen =[bond[0] for bond in self.reactant_bonds 
                                if bond[0] in self.fixed_atom or bond[1] in self.fixed_atom if self.atoms[bond[0]] == 1 and self.atoms[bond[1]] == 8]
            self.active_site_hydrogen =[bond[1] for bond in self.reactant_bonds 
                                if bond[0] in self.fixed_atom or bond[1] in self.fixed_atom if self.atoms[bond[1]] == 1 and self.atoms[bond[0]] == 8]
            self.active_site_hydrogen.extend(active_site_hydrogen)
        else:
            self.reactant_carbon = [idx for idx, reactant_atoms in enumerate(self.atoms) if reactant_atoms == 6]

    def gen_bond_can_form_and_bond_can_break(self, bonds_form_all):
        """
        Generate bond can form and bond can break for lewis acid or bronsted acid zeolite
        """

        bond_can_break = []
        if self.fixed_atom:
            if self.catalyst in ['SnBEA', 'WBEA', 'MoBEA']:
                bond_can_form = [bonds for bonds in bonds_form_all if bonds[0] not in self.fixed_atom and bonds[1] not in self.fixed_atom]
                # index start from 0
                # HR (reactant hydrogen, -OH) <---> OA (active site oxygen)
                # Like active site recovery
                rh_ao = [(hydrogen, oxygen, 1) 
                        for hydrogen in self.reactant_hydrogen 
                        for oxygen in self.active_site_oxygen]

                # Lewis acid Sn with oxygen
                # Reactant oxygen (OR)<---> Sn
                ro_metal = [(oxygen, metal, 1)
                            for oxygen in self.reactant_oxygen 
                            for metal in self.active_site_metal]

                # HR (reactant hydrogen, -CH) <---> OA (active site oxygen)
                ch_ao = [(hydrogen, oxygen, 1)
                        for hydrogen in self.reactant_hydrogen_2
                        for oxygen in self.active_site_oxygen]

                # Add to bond can form
                bond_can_form.extend(rh_ao)
                bond_can_form.extend(ro_metal)
                bond_can_form.extend(ch_ao)

                # Copy the bond can form
                bond_can_form_copy = bond_can_form[:]

                # Remove the cracking and dehydrogenation
                possible_bond_1 = [(carbon_and_hydrogen, hydrogen, 1) 
                                    for carbon_and_hydrogen in self.reactant_carbon + self.reactant_hydrogen_2 
                                    for hydrogen in self.active_site_hydrogen]
                possible_bond = [(hydrogen, carbon_and_hydrogen, 1) 
                                for carbon_and_hydrogen in self.reactant_carbon + self.reactant_hydrogen_2 
                                for hydrogen in self.active_site_hydrogen]
                possible_bond.extend(possible_bond_1)

                for bond in bond_can_form_copy:
                    # If active site have two hydrogen, preventing them become hydrogen
                    if bond[0] in self.active_site_hydrogen and bond[1] in self.active_site_hydrogen:
                        bond_can_form.remove(bond)
                    # reactant -OH --> active site hydrogen
                    if (bond[0] in self.reactant_hydrogen and bond[1] in self.active_site_hydrogen) or (bond[0] in self.active_site_hydrogen and bond[1] in self.reactant_hydrogen):
                        bond_can_form.remove(bond)
                    if bond in possible_bond:
                        bond_can_form.remove(bond)

                for bonds in self.reactant_bonds:
                    # Remove the original C-H bond (prevent CH bond type greater than 1 C=H)
                    if (self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 1) or (self.atoms[bonds[1]] == 6 and self.atoms[bonds[0]] == 1):
                        bond_can_form.remove(bonds)
                        # Prevent the dehydrogenation from the same CH. For example CH3  --> CH + H2
                        if self.atoms[bonds[0]] == 6:
                            hydrogen = [i[1] for i in self.reactant_bonds if i[0] == bonds[0] if self.atoms[i[1]] == 1]
                            combs = combinations(hydrogen, 2)
                            for comb in combs:
                                if (comb[0], comb[1], 1) in bond_can_form:
                                    bond_can_form.remove((comb[0], comb[1], 1))
                        elif self.atoms[bonds[1]] == 6:
                            hydrogen = [i[1] for i in self.reactant_bonds if i[0] == bonds[0] if self.atoms[i[1]] == 1]
                            combs = combinations(hydrogen, 2)
                            for comb in combs:
                                if (comb[0], comb[1], 1) in bond_can_form:
                                    bond_can_form.remove((comb[0], comb[1], 1))

                    # Remove the C-O triple bond
                    elif (self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 8) or (self.atoms[bonds[1]] == 8 and self.atoms[bonds[0]] == 6):
                        if bonds[2] == 2:
                            bond_can_form.remove((bonds[0], bonds[1], 1))

                    # Remove the original O-H bond (prevent OH bond type greater than 1 O=H)
                    elif (self.atoms[bonds[0]] == 1 and self.atoms[bonds[1]] == 8) or (self.atoms[bonds[1]] == 8 and self.atoms[bonds[0]] == 1):
                        if (bonds[0], bonds[1], 1) in bond_can_form:
                            bond_can_form.remove((bonds[0], bonds[1], 1))
                        elif (bonds[1], bonds[0], 1) in bond_can_form:
                            bond_can_form.remove((bonds[1], bonds[0], 1))

                    # Remove the C-C quadra bond
                    elif self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 6:
                        if bonds[2] == 3:
                            bond_can_form.remove((bonds[0], bonds[1], 1))

                    # Create bond can form list
                    if bonds[0] not in self.fixed_atom or bonds[1] not in self.fixed_atom:
                        order = bonds[2]
                        while order > 1:
                            bond_can_break.append((bonds[0], bonds[1], order - 1))
                            order -= 1
                        bond_can_break.append(bonds)

                bond_can_form = sorted(bond_can_form, key = itemgetter(0))
                bond_can_break = sorted(bond_can_break, key = itemgetter(0))
            elif self.catalyst in ['HZSM5']:
                bond_can_form = [bonds for bonds in bonds_form_all if bonds[0] not in self.fixed_atom and bonds[1] not in self.fixed_atom]
                # index start from 0
                # HR (reactant hydrogen, -OH) <---> OA (active site oxygen)
                # Like active site recovery
                rh_ao = [(hydrogen, oxygen, 1) 
                        for hydrogen in self.reactant_hydrogen 
                        for oxygen in self.active_site_oxygen]

                # HR (reactant hydrogen, -CH) <---> OA (active site oxygen)
                ch_ao = [(hydrogen, oxygen, 1)
                        for hydrogen in self.reactant_hydrogen_2
                        for oxygen in self.active_site_oxygen]

                # Add to bond can form
                bond_can_form.extend(rh_ao)
                bond_can_form.extend(ch_ao)

                # Copy the bond can form
                bond_can_form_copy = bond_can_form[:]

                for bond in bond_can_form_copy:
                    # If active site have two hydrogen, preventing them become hydrogen
                    if bond[0] in self.active_site_hydrogen and bond[1] in self.active_site_hydrogen:
                        bond_can_form.remove(bond)

                for bonds in self.reactant_bonds:
                    # Remove the original C-H bond (prevent CH bond type greater than 1 C=H)
                    if (self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 1) or (self.atoms[bonds[1]] == 6 and self.atoms[bonds[0]] == 1):
                        bond_can_form.remove(bonds)
                        # Prevent the dehydrogenation from the same CH. For example CH3  --> CH + H2
                        if self.atoms[bonds[0]] == 6:
                            hydrogen = [i[1] for i in self.reactant_bonds if i[0] == bonds[0] if self.atoms[i[1]] == 1]
                            combs = combinations(hydrogen, 2)
                            for comb in combs:
                                if (comb[0], comb[1], 1) in bond_can_form:
                                    bond_can_form.remove((comb[0], comb[1], 1))
                        elif self.atoms[bonds[1]] == 6:
                            hydrogen = [i[1] for i in self.reactant_bonds if i[0] == bonds[0] if self.atoms[i[1]] == 1]
                            combs = combinations(hydrogen, 2)
                            for comb in combs:
                                if (comb[0], comb[1], 1) in bond_can_form:
                                    bond_can_form.remove((comb[0], comb[1], 1))

                    # Remove the C-O triple bond
                    elif (self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 8) or (self.atoms[bonds[1]] == 8 and self.atoms[bonds[0]] == 6):
                        if bonds[2] == 2:
                            bond_can_form.remove((bonds[0], bonds[1], 1))

                    # Remove the original O-H bond (prevent OH bond type greater than 1 O=H)
                    elif (self.atoms[bonds[0]] == 1 and self.atoms[bonds[1]] == 8) or (self.atoms[bonds[1]] == 8 and self.atoms[bonds[0]] == 1):
                        if (bonds[0], bonds[1], 1) in bond_can_form:
                            bond_can_form.remove((bonds[0], bonds[1], 1))
                        elif (bonds[1], bonds[0], 1) in bond_can_form:
                            bond_can_form.remove((bonds[1], bonds[0], 1))

                    # Remove the C-C quadra bond
                    elif self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 6:
                        if bonds[2] == 3:
                            bond_can_form.remove((bonds[0], bonds[1], 1))

                    # Create bond can form list
                    if bonds[0] not in self.fixed_atom or bonds[1] not in self.fixed_atom:
                        order = bonds[2]
                        while order > 1:
                            bond_can_break.append((bonds[0], bonds[1], order - 1))
                            order -= 1
                        bond_can_break.append(bonds)

                bond_can_form = sorted(bond_can_form, key = itemgetter(0))
                bond_can_break = sorted(bond_can_break, key = itemgetter(0))
        else:
            bond_can_form = bonds_form_all[:]
            # Remove the original C-H bond
            for bonds in self.reactant_bonds:
                order = bonds[2]
                if (self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 1) or (self.atoms[bonds[1]] == 6 and self.atoms[bonds[0]] == 1):
                    bond_can_form.remove(bonds)
                # Remove the C-O triple bond
                elif (self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 8) or (self.atoms[bonds[1]] == 8 and self.atoms[bonds[0]] == 6):
                    if bonds[2] > 2:
                        bond_can_form.remove((bonds[0], bonds[1], 1))
                # Remove the C-C quadra bond
                elif self.atoms[bonds[0]] == 6 and self.atoms[bonds[1]] == 6:
                    if bonds[2] > 3:
                        bond_can_form.remove((bonds[0], bonds[1], 1))
                while order > 1:
                    bond_can_break.append((bonds[0], bonds[1], order - 1))
                    order -= 1
                bond_can_break.append(bonds)
            bond_can_form = sorted(bond_can_form, key = itemgetter(0))
            bond_can_break = sorted(bond_can_break, key = itemgetter(0))

        return bond_can_form, bond_can_break

    def generateProducts(self):
        """
        Generate all possible products from the reactant under the constraints
        of breaking a maximum of `nbreak` and forming a maximum of `nform`
        bonds.
        """
        if self.nbreak > 3 or self.nform > 3:
            raise Exception('Breaking/forming bonds is limited to a maximum of 3')

        # Extract valences as a mutable sequence
        reactant_valences = [atom.OBAtom.GetExplicitValence() for atom in self.reac_mol]
        # Initialize set for storing bonds of products
        # A set is used to ensure that no duplicate products are added
        products_bonds = set()

        # Generate all possibilities for forming bonds
        natoms = len(self.atoms)
        bonds_form_all = [(atom1_idx, atom2_idx, 1)
                          for atom1_idx in range(natoms - 1)
                          for atom2_idx in range(atom1_idx + 1, natoms)]

        # Generate bond can form and bond can break
        bond_can_form, bond_can_break = self.gen_bond_can_form_and_bond_can_break(bonds_form_all)

        # Generate products
        bf_combinations = ((0, 1), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3))

        for bf in bf_combinations:
            if bf[0] <= self.nbreak and bf[1] <= self.nform:
                self._generateProductsHelper(
                    bf[0],
                    bf[1],
                    products_bonds,
                    self.reactant_bonds,
                    reactant_valences,
                    bond_can_form,
                    bond_can_break
                )

        if products_bonds:
            for bonds in products_bonds:
                # for SSM calculation
                break_bonds = [i for i in list(set(bonds) ^ set(self.reactant_bonds)) if i not in bonds]
                form_bonds = [i for i in list(set(bonds) ^ set(self.reactant_bonds)) if i in bonds]
                
                break_bonds_copy = break_bonds[:]
                form_bonds_copy = form_bonds[:]
                # deal with double bonds
                for i in form_bonds_copy:
                    if i[2] > 1 and i in bonds:
                        for j in break_bonds_copy:
                            if i[0] == j[0] and i[1] == j[1]:
                                break_bonds.remove(j)

                # deal with double bonds
                for i in break_bonds_copy:
                    if i[2] > 1 and i in self.reactant_bonds:
                        for j in form_bonds_copy:
                            if i[0] == j[0] and i[1] == j[1]:
                                form_bonds.remove(j)

                if self.check_bond_type(bonds):
                    mol = gen3D.makeMolFromAtomsAndBonds(self.atoms, bonds, spin=self.reac_mol.spin)
                    mol.setCoordsFromMol(self.reac_mol)
                    if self.check_bond_dissociation_energy_and_isomorphic_and_rings(bonds, break_bonds):
                        if self.use_inchi_key and mol.write('inchiKey').strip() not in self.reactant_inchikey:
                            # Remove the double bond in forming or breaking bond.
                            # Because in SSM or GSM double bond is only a little distance change.
                            # That is the double bond can't be a driving coordinate, ssm will automatically deal with this little distance change.
                            for i in form_bonds_copy:
                                for j in self.reactant_bonds:
                                    if i[0] == j[0] and i[1] == j[1]:
                                        if i in form_bonds:
                                            form_bonds.remove(i)

                            for i in break_bonds_copy:
                                for j in bonds:
                                    if i[0] == j[0] and i[1] == j[1]:
                                        if i in break_bonds:
                                            break_bonds.remove(i)
                            
                            self.add_bonds.append(form_bonds)
                            self.break_bonds.append(break_bonds)
                            self.prod_mols.append(mol)

    def check_bond_type(self, bonds):
        bond_type = {}
        for i in range(len(self.atoms)):
            num = 0
            for j in bonds:
                if j[0] == i or j[1] == i:
                    num += j[2]
            bond_type[i] = num
        if 0 in bond_type.values():
            return False
        else:
            for idx, i in enumerate(self.atoms):
                # For catalyst (SnBEA)
                if self.fixed_atom:
                    if idx not in self.fixed_atom:
                        if i == 6 and bond_type[idx] != 4:
                            return False
                        elif i == 8 and bond_type[idx] < 2:
                            return False
                        elif i == 8 and bond_type[idx] > 2:
                            if (idx, self.active_site_metal[0], 1) not in bonds and (self.active_site_metal[0], idx, 1) not in bonds:
                                return False
                        elif i == 14 and bond_type[idx] != 4:
                            return False
                    else:
                        # While the bronsted acid already have proton on the active site, then aborted.
                        if i == 8 and bond_type[idx] > 3:
                            return False
                # Normal case (organic reaction)
                else:
                    if i == 6 and bond_type[idx] != 4:
                        return False
                    elif i == 8 and bond_type[idx] != 2:
                        return False
                    elif i == 14 and bond_type[idx] != 4:
                        return False
                    elif i == 1 and bond_type[idx] != 1:
                        return False
            return True

    def check_bond_dissociation_energy_and_isomorphic_and_rings(self, bond_list, bbond_list):
        energy = 0.0

        reactant_graph = self.reac_mol_graph
        atoms = [Atom(atomic_symbol=reactant_graph.atoms[i].label) for i in range(reactant_graph.n_atoms)]
        product = Species(atoms)
        graph = nx.Graph()

        for i in range(reactant_graph.n_atoms):
            graph.add_node(i, atom_label=reactant_graph.atoms[i].label, stereo=False)

        if bond_list is not None:
            [graph.add_edge(bond[0], bond[1], pi=False, active=False) for bond in bond_list]
            product.graph = graph

        # Filter the isomorphic
        if is_isomorphic(reactant_graph.graph, product.graph):
            return False
        elif self.check_four_and_three_membered_rings(product.graph):
            return False
        else:
            # Filter the bond dissociation energy
            num = sum([bb[2] for bb in bbond_list])
            for break_bond in bbond_list:
                first_atom = reactant_graph.atoms[break_bond[0]].label
                second_atom = reactant_graph.atoms[break_bond[1]].label
                bond_type = break_bond[2]
                supported_element = ['C', 'N', 'H', 'O', 'S', 'Cl', 'Si']
                if first_atom not in supported_element or second_atom not in supported_element:
                    # use 100 instead
                    energy += 100
                else:
                    # consider if break double bond (2-->1 not break 2) then the bond dissociation use double bond energy - single bond energy
                    if num >= 3 and bond_type >= 2:
                        try:
                            energy += props.bond_dissociation_energy[first_atom, second_atom, bond_type]
                            energy -= props.bond_dissociation_energy[first_atom, second_atom, bond_type - 1]
                        except:
                            energy += props.bond_dissociation_energy[second_atom, first_atom, bond_type]
                            energy -= props.bond_dissociation_energy[second_atom, first_atom, bond_type - 1]
                    else:
                        try:
                            energy += props.bond_dissociation_energy[first_atom, second_atom, bond_type]
                        except:
                            energy += props.bond_dissociation_energy[second_atom, first_atom, bond_type]
                            
            if energy / _constants.CAL2J > self.bond_dissociation_cutoff:
                return False
            else:
                return True

    def check_four_and_three_membered_rings(self, product):
        # Filter the 3&4 membered ring
        rings = nx.cycle_basis(product)
        for ring in rings:
            if len(ring) < 5 and set(self.reactant_carbon) >= set(ring):
                return True
        return False

    @property
    def get_prods(self):
        return self.prod_mols

    @property
    def get_add_bonds(self):
        return self.add_bonds

    @property
    def get_break_bonds(self):
        return self.break_bonds

    def _generateProductsHelper(self, nbreak, nform, products, bonds, valences, bond_can_form, bond_can_break, bonds_broken=None):
        """
        Generate products recursively given the number of bonds that should be
        broken and formed, a set for storing the products, a sequence of atoms,
        of bonds, and of valences. `bond_can_form` should contain a tuple of
        tuples of bonds that contains all possibilities for forming bonds.
        Nothing is returned, but formed products are added to `products`.
        """
        if bonds_broken is None:
            bonds_broken = []

        if nbreak == 0 and nform == 0:
            products.add((tuple(sorted(bonds))))
            # if all(bonds_broken[num][0] not in self.fixed_atom or bonds_broken[num][1] not in self.fixed_atom for num in range(len(bonds_broken))):
            #     products.add((tuple(sorted(bonds))))

        if nbreak > 0:
            # Break bond
            for bond_break in bond_can_break:
                if bond_break not in bonds:
                    continue
                bond_break_idx = bonds.index(bond_break)
                valences_break = self.changeValences(valences, bond_break, -1)
                bonds_break = self.breakBond(bonds, bond_break_idx)

                # Keep track of bonds that have been broken
                if bond_break_idx == 0:
                    bonds_broken.append(bond_break)
                else:
                    bonds_broken[-1] = bond_break

                # Call function recursively to break next bond
                self._generateProductsHelper(
                    nbreak - 1,
                    nform,
                    products,
                    bonds_break,
                    valences_break,
                    bond_can_form,
                    bond_can_break,
                    bonds_broken
                )
            # Remove last bond that has been broken after loop terminates
            del bonds_broken[-1]
        elif nform > 0:
            # Form bond
            for bond_form in bond_can_form:
                # Do not add bond if it has previously been broken
                if bond_form[:2] in [bond[:2] for bond in bonds_broken]:
                    continue
                # Form new bond and catch exception if it violates constraints
                try:
                    valences_form = self.changeValences(valences, bond_form, 1)
                    bonds_form = self.formBond(bonds, bond_form)
                except StructureError:
                    continue

                # Call function recursively to form next bond
                self._generateProductsHelper(
                    nbreak,
                    nform - 1,
                    products,
                    bonds_form,
                    valences_form,
                    bond_can_form,
                    bond_can_break,
                    bonds_broken
                )

    @staticmethod
    def breakBond(bonds, break_idx):
        """
        Break a bond given a tuple of tuples of bonds in the form:
            (atom index 1, atom index 2, bond type)
        The bond at index `break_idx` is broken. A tuple of tuples of updated
        bonds is returned.
        """
        # Break double or triple bond
        if bonds[break_idx][2] > 1:
            return bonds[:break_idx] + (bonds[break_idx][:2] + (bonds[break_idx][2] - 1,),) + bonds[(break_idx + 1):]

        # Break single bond
        return bonds[:break_idx] + bonds[(break_idx + 1):]

    @staticmethod
    def formBond(bonds, new_bond):
        """
        Form a bond given a tuple of tuples of bonds in the form:
            (atom index 1, atom index 2, bond type)
        If a bond already exists in `bonds`, the bond order is incremented. If
        a bond between atoms that are not yet connected is to be formed, then a
        new bond is added to the tuple of bonds. The bond to be added is
        specified in `new_bond`. A tuple of tuples of updated bonds is
        returned.
        """
        # Ensure that only one bond is added at a time
        assert new_bond[2] == 1

        try:  # Check if bond exists as single bond
            idx = bonds.index(new_bond)
        except ValueError:
            try:  # Check if bond exists as double bond
                idx = bonds.index(new_bond[:2] + (2,))
            except ValueError:
                try:  # Check if bond exists as triple bond
                    idx = bonds.index(new_bond[:2] + (3,))
                except ValueError:  # Add new bond if it does not exist yet
                    return bonds + (new_bond,)
                else:  # Raise exception if trying to exceed triple bond
                    raise StructureError('Bond type cannot be higher than triple bond for bond {}'.format(bonds[idx]))
            else:  # Return bonds with double bond increased to triple bond
                return bonds[:idx] + (bonds[idx][:2] + (3,),) + bonds[(idx + 1):]
        else:  # Return bonds with single bond increased to double bond
            return bonds[:idx] + (bonds[idx][:2] + (2,),) + bonds[(idx + 1):]

    def changeValences(self, valences, bond, inc):
        """
        Update the valences corresponding to each atom in `self.atoms` given
        the current valences in `valences`, the bond that is being affected,
        and the increment (typically 1 or -1). The maximum valences for each
        atom type are respected and an error is raised if they would be
        violated. A sequence of updated valences is returned.
        """
        # Create copy of current valences
        valences_temp = valences[:]

        # Check for invalid operation
        if inc < 0 and (valences_temp[bond[0]] < inc or valences_temp[bond[1]] < inc):
            raise Exception('Cannot decrease valence below zero-valence')

        # Change valences of both atoms participating in bond
        valences_temp[bond[0]] += inc
        valences_temp[bond[1]] += inc

        # Check if maximum valences are exceeded
        element0 = ELEMENT_TABLE.from_atomic_number(self.atoms[bond[0]])
        element1 = ELEMENT_TABLE.from_atomic_number(self.atoms[bond[1]])
        if valences_temp[bond[0]] > element0.max_bonds:
            raise StructureError('Maximum valence on atom {} exceeded'.format(bond[0]))
        if valences_temp[bond[1]] > element1.max_bonds:
            raise StructureError('Maximum valence on atom {} exceeded'.format(bond[1]))

        # Return valid valences
        return valences_temp
