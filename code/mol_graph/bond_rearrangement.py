# standard library imports
import os

# third party
import itertools

# local application imports
from atoms import get_maximal_valance
from geom import get_neighbour_list
from geom import get_points_on_sphere
from graph import get_bond_type_list
from graph import get_fbonds
from graph import is_isomorphic
from graph import connected_components


def get_bond_rearrangs(reactant, product):
    """For a reactant and product (complex) find the set of breaking and
    forming bonds that will turn reactants into products. This works by
    determining the types of bonds that have been made/broken (i.e CH) and
    then only considering rearrangements involving those bonds.

    Arguments:
        reactant (autode.complex.ReactantComplex):
        product (autode.complex.ProductComplex):
        name (str): Name of the reaction

    Returns:
        list: list of bond rearrang objects linking reacs and prods
    """
    """
    if os.path.exists(f'{name}_bond_rearrangs.txt'):
        return get_bond_rearrangs_from_file(f'{name}_bond_rearrangs.txt')

    if is_isomorphic(reactant.graph, product.graph) and product.n_atoms > 3:
        return None
    """
    possible_bond_rearrs = []

    reac_bond_dict = get_bond_type_list(reactant.graph)
    prod_bond_dict = get_bond_type_list(product.graph)
    # list of bonds where this type of bond (e.g C-H) has less bonds in
    # products than reactants
    all_possible_bbonds = []

    # list of bonds that can be formed of this bond type. This is only used
    # if there is only one type of bbond, so can be overwritten for each new
    # type of bbond
    bbond_atom_type_fbonds = None

    # list of bonds where this type of bond (e.g C-H) has more bonds in
    #  products than reactants
    all_possible_fbonds = []

    # list of bonds that can be broken of this bond type. This is only used
    # if there is only one type of fbond, so can be overwritten for each new
    # type of fbond
    fbond_atom_type_bbonds = None

    # list of bonds where this type of bond (e.g C-H) has the same number of
    # bonds in products and reactants
    possible_bbond_and_fbonds = []

    for reac_key, reac_bonds in reac_bond_dict.items():
        prod_bonds = prod_bond_dict[reac_key]
        possible_fbonds = get_fbonds(reactant.graph, reac_key)
        if len(prod_bonds) < len(reac_bonds):
            all_possible_bbonds.append(reac_bonds)
            bbond_atom_type_fbonds = possible_fbonds
        elif len(prod_bonds) > len(reac_bonds):
            all_possible_fbonds.append(possible_fbonds)
            fbond_atom_type_bbonds = reac_bonds
        else:
            if len(reac_bonds) != 0:
                possible_bbond_and_fbonds.append([reac_bonds, possible_fbonds])

    # The change in the number of bonds is > 0 as in the reaction
    # initialisation reacs/prods are swapped if this is < 0
    delta_n_bonds = reactant.graph.number_of_edges() - product.graph.number_of_edges()
    if delta_n_bonds == 0:
        funcs = [get_fbonds_bbonds_1b1f, get_fbonds_bbonds_2b2f]
    elif delta_n_bonds == 1:
        funcs = [get_fbonds_bbonds_1b, get_fbonds_bbonds_2b1f]
    elif delta_n_bonds == 2:
        funcs = [get_fbonds_bbonds_2b]
    else:
        return None

    for func in funcs:
        possible_bond_rearrs = func(reactant, product, possible_bond_rearrs,
                                    all_possible_bbonds,
                                    all_possible_fbonds,
                                    possible_bbond_and_fbonds,
                                    bbond_atom_type_fbonds,
                                    fbond_atom_type_bbonds)

        if len(possible_bond_rearrs) > 0:
            # This function will return with the first bond rearrangement
            # that leads to products
            n_bond_rearrangs = len(possible_bond_rearrs)
            if n_bond_rearrangs > 1:
                possible_bond_rearrs = strip_equivalent_bond_rearrangs(
                    reactant, possible_bond_rearrs)

            # save_bond_rearrangs_to_file(possible_bond_rearrs,
                # filename=f'{name}_bond_rearrangs.txt')

            return possible_bond_rearrs

    return None


def save_bond_rearrangs_to_file(bond_rearrangs, filename='bond_rearrangs.txt'):
    with open(filename, 'w') as file:
        for bond_rearrang in bond_rearrangs:
            print('fbonds', file=file)
            for fbond in bond_rearrang.fbonds:
                print(*fbond, file=file)
            print('bbonds', file=file)
            for bbond in bond_rearrang.bbonds:
                print(*bbond, file=file)
            print('end', file=file)

    return None


def get_bond_rearrangs_from_file(filename='bond_rearrangs.txt'):

    if not os.path.exists(filename):
        return None

    bond_rearrangs = []

    with open(filename, 'r') as br_file:
        fbonds_block = False
        bbonds_block = True
        fbonds = []
        bbonds = []
        for line in br_file:
            if 'fbonds' in line:
                fbonds_block = True
                bbonds_block = False
            if 'bbonds' in line:
                fbonds_block = False
                bbonds_block = True
            if fbonds_block and len(line.split()) == 2:
                atom_id_string = line.split()
                fbonds.append((int(atom_id_string[0]), int(atom_id_string[1])))
            if bbonds_block and len(line.split()) == 2:
                atom_id_string = line.split()
                bbonds.append((int(atom_id_string[0]), int(atom_id_string[1])))
            if 'end' in line:
                bond_rearrangs.append(BondRearrangement(forming_bonds=fbonds,
                                                        breaking_bonds=bbonds))
                fbonds = []
                bbonds = []

    return bond_rearrangs


def generate_rearranged_graph(graph, fbonds, bbonds):
    """Generate a rearranged graph by breaking bonds (edge) and forming others (edge)

    Arguments:
        graph (nx.Graph): reactant graph
        fbonds (list(tuple)): list of bonds to be made
        bbonds (list(tuple)): list of bonds to be broken

    Returns:
        nx.Graph: rearranged graph
    """
    rearranged_graph = graph.copy()
    for fbond in fbonds:
        rearranged_graph.add_edge(*fbond)
    for bbond in bbonds:
        rearranged_graph.remove_edge(*bbond)

    return rearranged_graph


def add_bond_rearrangment(bond_rearrangs, reactant, product, fbonds, bbonds):
    """For a possible bond rearrangement, sees if the products are made, and
    adds it to the bond rearrang list if it does

    Arguments:
        bond_rearrangs (list): list of working bond rearrangments
        reactant (molecule object): reactant complex
        product (molecule object): product complex
        fbonds (list of tuples): list of bonds to be made
        bbonds (list of tuples): list of bonds to be broken

    Returns:
        list: updated list of working bond rearrangments
    """

    # Check that the bond rearrangement doesn't exceed standard atom valances
    bbond_atoms = [atom for bbond in bbonds for atom in bbond]
    for fbond in fbonds:
        for atom in fbond:
            atom_label = reactant.atoms[atom].label
            if reactant.graph.degree(atom) == get_maximal_valance(atom_label) and atom not in bbond_atoms:
                # If we are here then there is at least one atom that will
                # exceed it's maximal valance, therefore
                # we don't need to run isomorphism
                return bond_rearrangs

    rearranged_graph = generate_rearranged_graph(
        reactant.graph, fbonds=fbonds, bbonds=bbonds)
    if is_isomorphic(rearranged_graph, product.graph):
        ordered_fbonds = []
        ordered_bbonds = []
        for fbond in fbonds:
            if fbond[0] < fbond[1]:
                ordered_fbonds.append((fbond[0], fbond[1]))
            else:
                ordered_fbonds.append((fbond[1], fbond[0]))
        for bbond in bbonds:
            if bbond[0] < bbond[1]:
                ordered_bbonds.append((bbond[0], bbond[1]))
            else:
                ordered_bbonds.append((bbond[1], bbond[0]))
        ordered_fbonds.sort()
        ordered_bbonds.sort()
        bond_rearrangs.append(BondRearrangement(forming_bonds=ordered_fbonds,
                                                breaking_bonds=ordered_bbonds))

    return bond_rearrangs


def get_fbonds_bbonds_1b(reactant, product, possible_bond_rearrangs, all_possible_bbonds, all_possible_fbonds, possible_bbond_and_fbonds, bbond_atom_type_fbonds, fbond_atom_type_bbonds):

    for bbond in all_possible_bbonds[0]:
        # break one bond
        possible_bond_rearrangs = add_bond_rearrangment(
            possible_bond_rearrangs, reactant, product, fbonds=[], bbonds=[bbond])

    return possible_bond_rearrangs


def get_fbonds_bbonds_2b(reactant, product, possible_bond_rearrangs, all_possible_bbonds, all_possible_fbonds, possible_bbond_and_fbonds, bbond_atom_type_fbonds, fbond_atom_type_bbonds):

    if len(all_possible_bbonds) == 1:
        # break two bonds of the same type
        for bbond1, bbond2 in itertools.combinations(all_possible_bbonds[0], 2):
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[], bbonds=[bbond1, bbond2])

    elif len(all_possible_bbonds) == 2:
        # break two bonds of different types
        for bbond1, bbond2 in itertools.product(all_possible_bbonds[0], all_possible_bbonds[1]):
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[], bbonds=[bbond1, bbond2])

    return possible_bond_rearrangs


def get_fbonds_bbonds_1b1f(reactant, product, possible_bond_rearrangs, all_possible_bbonds, all_possible_fbonds, possible_bbond_and_fbonds, bbond_atom_type_fbonds, fbond_atom_type_bbonds):

    if len(all_possible_bbonds) == 1 and len(all_possible_fbonds) == 1:
        # make and break a bond of different types
        for fbond, bbond in itertools.product(all_possible_fbonds[0], all_possible_bbonds[0]):
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond], bbonds=[bbond])

    elif len(all_possible_bbonds) == 0 and len(all_possible_fbonds) == 0:
        # make and break a bond of the same type
        for bbonds, fbonds in possible_bbond_and_fbonds:
            for bbond, fbond in itertools.product(bbonds, fbonds):
                possible_bond_rearrangs = add_bond_rearrangment(
                    possible_bond_rearrangs, reactant, product, fbonds=[fbond], bbonds=[bbond])

    return possible_bond_rearrangs


def get_fbonds_bbonds_2b1f(reactant, product, possible_bond_rearrangs, all_possible_bbonds, all_possible_fbonds, possible_bbond_and_fbonds, bbond_atom_type_fbonds, fbond_atom_type_bbonds):

    if len(all_possible_bbonds) == 2 and len(all_possible_fbonds) == 1:
        # make a bond and break two bonds, all of different types
        for fbond, bbond1, bbond2 in itertools.product(all_possible_fbonds[0], all_possible_bbonds[0], all_possible_bbonds[1]):
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond], bbonds=[bbond1, bbond2])

    elif len(all_possible_bbonds) == 1 and len(all_possible_fbonds) == 1:
        # make a bond of one type, break two bonds of another type
        for fbond, (bbond1, bbond2) in itertools.product(all_possible_fbonds[0], itertools.combinations(all_possible_bbonds[0], 2)):
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond], bbonds=[bbond1, bbond2])

    elif len(all_possible_bbonds) == 1 and len(all_possible_fbonds) == 0:
        for bbonds, fbonds in possible_bbond_and_fbonds:
            # make and break a bond of one type, break a bond of a different type
            for fbond, bbond1, bbond2 in itertools.product(fbonds, all_possible_bbonds[0], bbonds):
                possible_bond_rearrangs = add_bond_rearrangment(
                    possible_bond_rearrangs, reactant, product, fbonds=[fbond], bbonds=[bbond1, bbond2])

        for fbond, (bbond1, bbond2) in itertools.product(bbond_atom_type_fbonds, itertools.combinations(all_possible_bbonds[0], 2)):
            # make and break two bonds, all of the same type
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond], bbonds=[bbond1, bbond2])

    return possible_bond_rearrangs


def get_fbonds_bbonds_2b2f(reactant, product, possible_bond_rearrangs, all_possible_bbonds, all_possible_fbonds, possible_bbond_and_fbonds, bbond_atom_type_fbonds, fbond_atom_type_bbonds):

    if len(all_possible_bbonds) == 2 and len(all_possible_fbonds) == 2:
        # make two bonds and break two bonds, all of different types
        for fbond1, fbond2, bbond1, bbond2 in itertools.product(all_possible_fbonds[0], all_possible_fbonds[1], all_possible_bbonds[0], all_possible_bbonds[1]):
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

    elif len(all_possible_bbonds) == 2 and len(all_possible_fbonds) == 1:
        # make two bonds of the same type, break two bonds of different types
        for bbond1, bbond2, (fbond1, fbond2) in itertools.product(all_possible_bbonds[0], all_possible_bbonds[1], itertools.combinations(all_possible_fbonds[0], 2)):
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

    elif len(all_possible_bbonds) == 1 and len(all_possible_fbonds) == 2:
        # make two bonds of different types, break two bonds of the same type
        for fbond1, fbond2, (bbond1, bbond2) in itertools.product(all_possible_fbonds[0], all_possible_fbonds[1], itertools.combinations(all_possible_bbonds[0], 2)):
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

    elif len(all_possible_bbonds) == 1 and len(all_possible_fbonds) == 1:
        for (fbond1, fbond2), (bbond1, bbond2) in itertools.product(itertools.combinations(all_possible_fbonds[0], 2), itertools.combinations(all_possible_bbonds[0], 2)):
            # make two bonds of the same type, break two bonds of another type
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

        for bbonds, fbonds in possible_bbond_and_fbonds:
            # make one bonds of one type, break one bond of another type, make and break a bond of a third type
            for fbond1, fbond2, bbond1, bbond2 in itertools.product(all_possible_fbonds[0], fbonds, all_possible_bbonds[0], bbonds):
                possible_bond_rearrangs = add_bond_rearrangment(
                    possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

        for fbond1, fbond2, (bbond1, bbond2) in itertools.product(all_possible_fbonds[0], bbond_atom_type_fbonds, itertools.combinations(all_possible_bbonds[0], 2)):
            # make a bond of one type, make and break two bonds of another type
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

        for bbond1, bbond2, (fbond1, fbond2) in itertools.product(all_possible_bbonds[0], fbond_atom_type_bbonds, itertools.combinations(all_possible_fbonds[0], 2)):
            # break a bond of one type, make two and break one bond of another type
            possible_bond_rearrangs = add_bond_rearrangment(
                possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

    elif len(all_possible_bbonds) == 0 and len(all_possible_fbonds) == 0:
        for (bbonds1, fbonds1), (bbonds2, fbonds2) in itertools.combinations(possible_bbond_and_fbonds, 2):
            # make and break a bond of one type, make and break a bond of another type
            for fbond1, bbond1, fbond2, bbond2 in itertools.product(fbonds1, bbonds1, fbonds2, bbonds2):
                possible_bond_rearrangs = add_bond_rearrangment(
                    possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

        for bbonds, fbonds in possible_bbond_and_fbonds:
            # make two and break two bonds, all of the same type
            for (fbond1, fbond2), (bbond1, bbond2) in itertools.product(itertools.combinations(fbonds, 2), itertools.combinations(bbonds, 2)):
                possible_bond_rearrangs = add_bond_rearrangment(
                    possible_bond_rearrangs, reactant, product, fbonds=[fbond1, fbond2], bbonds=[bbond1, bbond2])

    return possible_bond_rearrangs


def strip_equivalent_bond_rearrangs(mol, possible_bond_rearrs, depth=6):
    """Remove any bond rearrangement from possible_bond_rearrs for which
    there is already an equivalent in the unique_bond_rearrangements list

    Arguments:
        mol (molecule object): reactant object
        possible_bond_rearrs (list(object)): list of BondRearrangement objects

    Keyword Arguments:
        depth (int): Depth of neighbour list that must be identical for a set
               of atoms to be considered equivalent (default: {6})

    Returns:
        (list(BondRearrangement)): stripped list of BondRearrangement objects
    """

    unique_bond_rearrs = []

    for bond_rearr in possible_bond_rearrs:
        bond_rearrang_is_unique = True

        # Compare bond_rearrang to all those already considered to be unique,
        for unique_bond_rearrang in unique_bond_rearrs:

            if (unique_bond_rearrang.get_active_atom_neighbour_lists(mol=mol, depth=depth) ==
                    bond_rearr.get_active_atom_neighbour_lists(mol=mol, depth=depth)):
                bond_rearrang_is_unique = False

        if bond_rearrang_is_unique:
            unique_bond_rearrs.append(bond_rearr)

    return unique_bond_rearrs


class BondRearrangement:

    def __str__(self):
        return '_'.join(f'{bond[0]}-{bond[1]}' for bond in self.all)

    def get_active_atom_neighbour_lists(self, mol, depth):
        """
        Get neighbour lists of all the active atoms in the molecule
        (reactant complex)

        Arguments:
            mol (autode.species.Species):
            depth (int): Depth of the neighbour list to consider

        Returns:
            (list(list(int))):
        """
        connected_molecules = connected_components(mol.graph)
        n_molecules = len(connected_molecules)

        def shift_molecules(vectors):
            for i, molecule_nodes in enumerate(connected_molecules):
                for j in molecule_nodes:
                    mol.atoms[j].translate(vec=vectors[i])

        # For every molecule in the complex shift so they are far away, thus
        # the neighbour lists only include atoms in the same molecule
        shift_vectors = [
            100 * vec for vec in get_points_on_sphere(n_points=n_molecules+1)]
        shift_molecules(vectors=shift_vectors)

        # Calculate the neighbour lists while the molecules are all far away
        if self.active_atom_nl is None:
            self.active_atom_nl = [get_neighbour_list(species=mol, atom_i=atom)[
                :depth] for atom in self.active_atoms]
        # Shift the molecules back to where they were
        shift_molecules(vectors=[-vector for vector in shift_vectors])

        return self.active_atom_nl

    def _set_active_atom_list(self, bonds, ls):

        for bond in bonds:
            for atom in bond:
                if atom not in ls:
                    ls.append(atom)
                if atom not in self.active_atoms:
                    self.active_atoms.append(atom)

        return None

    def __eq__(self, other):
        return self.fbonds == other.fbonds and self.bbonds == other.bbonds

    def __init__(self, forming_bonds=None, breaking_bonds=None):
        """
        Bond rearrangement

        Keyword Arguments:
            forming_bonds (list(tuple(int))): List of atom pairs that are
                        forming in this reaction

            breaking_bonds (list(tuple(int))): List of atom pairs that are
                           breaking in the reaction
        """

        self.fbonds = forming_bonds if forming_bonds is not None else []
        self.bbonds = breaking_bonds if breaking_bonds is not None else []

        self.n_fbonds = len(self.fbonds)
        self.n_bbonds = len(self.bbonds)

        self.active_atoms = []
        self.fatoms = []
        self.batoms = []
        self.active_atom_nl = None
        self.all = self.fbonds + self.bbonds

        self._set_active_atom_list(self.fbonds, self.fatoms)
        self._set_active_atom_list(self.bbonds, self.batoms)
