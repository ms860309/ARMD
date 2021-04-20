# standard library imports
import os
import signal

# third party
import networkx as nx
import numpy as np
from scipy.spatial import distance_matrix
import itertools
from networkx.algorithms import isomorphism
from copy import deepcopy


def make_graph(species, rel_tolerance=0.25, bond_list=None,
               allow_invalid_valancies=False):
    """
    Make the molecular graph from the 'bonds' determined on a distance criteria
    or a smiles parser object. All attributes default to false

    Nodes attributes:
        (0) atom_label: Atomic symbol of this atom
        (1) stereo: Is this atom part of some stereochemistry e.g. R/S or E/Z

    Edge attributes:
        (1) pi: Is this bond a pi bond. If it is then there should be no
                rotation the bond axis in conformer generation
        (2) active: Is this bond being made/broken
                   (applies only to TransitionState objects)

    Arguments:
        species (autode.species.Species):

    Keyword Arguments:
        rel_tolerance (float):
        bond_list (list(tuple)):
        allow_invalid_valancies (bool):
    """

    graph = nx.Graph()

    # Add the atoms to the graph all are initially assumed not to be
    # stereocenters
    for i in range(species.n_atoms):
        graph.add_node(i, atom_label=species.atoms[i].label, stereo=False)

    # If bonds are specified then add edges to the graph and return
    if bond_list is not None:
        [graph.add_edge(bond[0], bond[1], pi=False, active=False)
         for bond in bond_list]
        species.graph = graph
        return None

    else:
        # Loop over the unique pairs of atoms and add 'bonds'
        coordinates = species.get_coordinates()
        dist_mat = distance_matrix(coordinates, coordinates)

        for i in get_atom_ids_sorted_type(species):

            # Iterate through the closest atoms to atom i
            for j in np.argsort(dist_mat[i]):

                if i == j:
                    # Don't bond atoms to themselves
                    continue

                avg_bond_length = get_avg_bond_length(atom_i_label=species.atoms[i].label,
                                                      atom_j_label=species.atoms[j].label)

                # If the distance between atoms i and j are less or equal to 1.2x average length add a 'bond'
                if dist_mat[i, j] <= avg_bond_length * (1.0 + rel_tolerance) and (i, j) not in graph.edges:
                    graph.add_edge(i, j, pi=False, active=False)

    species.graph = graph
    set_graph_attributes(species)

    if not allow_invalid_valancies:
        remove_bonds_invalid_valancies(species)

    return None


def get_bond_type_list(graph):
    """Finds the types (i.e CH) of all the bonds in a molecular graph
    Arguments:
        graph (nx.Graph): Molecular graph
    Returns:
        bond_list_dict (dict): key = bond type, value = list of bonds of this type
    """
    bond_list_dict = {}
    atom_types = set()

    for _, atom_label in graph.nodes.data('atom_label'):
        atom_types.add(atom_label)

    ordered_atom_labels = sorted(atom_types)

    for index, atom_label in enumerate(ordered_atom_labels):
        for i in range(index, len(ordered_atom_labels)):
            key = atom_label + ordered_atom_labels[i]
            bond_list_dict[key] = []

    for bond in graph.edges:
        atom_i_label = graph.nodes[bond[0]]['atom_label']
        atom_j_label = graph.nodes[bond[1]]['atom_label']
        key1, key2 = atom_i_label + atom_j_label, atom_j_label + atom_i_label

        if key1 in bond_list_dict.keys():
            bond_list_dict[key1].append(bond)
        elif key2 in bond_list_dict.keys():
            bond_list_dict[key2].append(bond)

    return bond_list_dict


def get_fbonds(graph, key):
    """Get all the possible forming bonds of a certain type
    Arguments:
        graph (nx.Graph): graph object of a molecule
        key (str): string representing the bond type to be examined
    Returns:
        list: list of bonds that can be made of this type
    """
    possible_fbonds = []
    bonds = list(graph.edges)
    for atom_i in graph.nodes:
        for atom_j in graph.nodes:
            if atom_i < atom_j:
                if not (atom_i, atom_j) in bonds and not (atom_j, atom_i) in bonds:
                    bond = (atom_i, atom_j)
                    atom_i_label = graph.nodes[bond[0]]['atom_label']
                    atom_j_label = graph.nodes[bond[1]]['atom_label']
                    key1, key2 = atom_i_label + atom_j_label, atom_j_label + atom_i_label
                    if key1 == key or key2 == key:
                        possible_fbonds.append(bond)

    return possible_fbonds


def connected_components(graph):
    return list(nx.connected_components(graph))


def get_atom_ids_sorted_type(species):
    """
    Get a list of atom ids sorted by increasing atomic weight, useful for when
     a molecular graph depends on the order
    of atoms in what will be considered bonded

    Arguments:
        species (autode.species.Species):

    Returns:
        (list(int)):
    """
    return sorted(list(range(species.n_atoms)), key=lambda i: get_atomic_weight(atom_label=species.atoms[i].label))


def set_graph_attributes(species):
    """
    For a molecular species set the π bonds and stereocentres in the molecular graph.

    Arguments:
        species (autode.species.Species):
    """

    for bond in species.graph.edges:
        atom_i, atom_j = bond

        if all([is_pi_atom(atom_label=species.atoms[atom].label, valency=species.graph.degree[atom]) for atom in bond]):
            # TODO fix this for alternating single and double bonds, currently all shown as pi
            species.graph.edges[atom_i, atom_j]['pi'] = True

    # List of atom indexes that are rings in the species
    rings = find_cycles(species.graph)

    for (i, j) in species.graph.edges:

        if species.graph.edges[(i, j)]['pi'] is False:
            continue

        if any(i in ring for ring in rings):
            # The ring should define the stereochemistry of this pi bond
            continue

        if is_chiral_pi_bond(species, bond=(i, j)):
            species.graph.nodes[i]['stereo'] = True
            species.graph.nodes[j]['stereo'] = True

    for i in range(species.n_atoms):
        if is_chiral_atom(species, atom_index=i):
            species.graph.nodes[i]['stereo'] = True

    return None


def find_cycles(graph):
    """Finds all the cycles in a graph

    Arguments:
        graph (nx.Graph): the molecular graph

    Returns:
        list(list): each list has the atoms in a cycle
    """
    return nx.cycle_basis(graph)


def is_chiral_atom(species, atom_index):
    """Determine if an atom is chiral, by seeing if any of the bonded groups
    are the same"""
    neighbours = list(species.graph.neighbors(atom_index))

    if len(neighbours) != 4:
        return False

    graphs = []
    for neighbour in neighbours:
        graph = species.graph.copy()
        graph.remove_edge(atom_index, neighbour)
        split_subgraphs = get_separate_subgraphs(graph)
        graphs.append([subgraph for subgraph in split_subgraphs if neighbour in list(
            subgraph.nodes())][0])

    for graph1, graph2 in itertools.combinations(graphs, 2):
        if is_isomorphic(graph1, graph2, ignore_active_bonds=True):
            return False

    return True


def reac_graph_to_prod_graph(reac_graph, bond_rearrang):
    """Makes the graph of the product from the reactant and the bond rearrang,
    so it has the indices of the reactant
    Arguments:
        reac_graph (nx.Graph): graph of the reactant
        bond_rearrang (autode.bond_rearrangement.BondRearrangement): the bond
                       rearrang linking reacs and prods
    Returns:
        nx.Graph: graph of the product with each atom indexed as in the reactants
    """
    prod_graph = deepcopy(reac_graph)
    for fbond in bond_rearrang.fbonds:
        prod_graph.add_edge(*fbond)
    for bbond in bond_rearrang.bbonds:
        prod_graph.remove_edge(*bbond)
    return prod_graph


def get_separate_subgraphs(graph):
    """Find all the unconnected graphs in a graph

    Arguments:
        graph (nx.Graph): graph

    Returns:
        list: list of graphs separate graphs
    """
    return [graph.subgraph(c).copy() for c in nx.connected_components(graph)]


def is_isomorphic(graph1, graph2, ignore_active_bonds=False, timeout=5):
    """Check whether two NX graphs are isomorphic. Contains a timeout because
    the gm.is_isomorphic() method occasionally gets stuck

    Arguments:
        graph1 (nx.Graph): graph 1
        graph2 (nx.Graph): graph 2

    Keyword Arguments:
        ignore_active_bonds (bool):
        timeout (float): Timeout in seconds

    Returns:
        (bool): if the graphs are isomorphic
    """

    if ignore_active_bonds:
        graph1, graph2 = get_graphs_ignoring_active_edges(graph1, graph2)

    if not isomorphism.faster_could_be_isomorphic(graph1, graph2):
        return False

    # Always match on atom types
    node_match = isomorphism.categorical_node_match('atom_label', 'C')

    if ignore_active_bonds:
        gm = isomorphism.GraphMatcher(graph1, graph2,
                                      node_match=node_match)

    else:
        # Also match on edges
        edge_match = isomorphism.categorical_edge_match('active', False)
        gm = isomorphism.GraphMatcher(graph1, graph2,
                                      node_match=node_match,
                                      edge_match=edge_match)

    # NX can hang here for not very large graphs, so kill after a timeout

    def handler(signum, frame):
        raise TimeoutError

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(int(timeout))
    try:
        result = gm.is_isomorphic()
        # Cancel the timer
        signal.alarm(0)
        return result

    except TimeoutError:
        return False


def get_graphs_ignoring_active_edges(graph1, graph2):
    """
    Remove any active edges that are in either graph1 or graph2 from both graphs
    Arguments:
        graph1 (nx.Graph):
        graph2 (nx.Graph):

    Returns:
        (tuple(nx.Graph))
    """
    graph1_no_ae, graph2_no_ae = graph1.copy(), graph2.copy()

    # Iterate through the pairs removing any active edges from both ga and gb
    for (ga, gb) in [(graph1_no_ae, graph2_no_ae), (graph2_no_ae, graph1_no_ae)]:

        for (i, j) in [edge for edge in ga.edges if ga.edges[edge]['active'] is True]:
            ga.remove_edge(i, j)

            if (i, j) in gb.edges:
                gb.remove_edge(i, j)

    return graph1_no_ae, graph2_no_ae


def remove_bonds_invalid_valancies(species):
    """
    Remove invalid valencies for atoms that exceed their maximum valencies e.g.
    H should have no more than 1 'bond'

    Arguments:
        species (autode.species.Species):
    """

    for i in species.graph.nodes:

        max_valance = get_maximal_valance(atom_label=species.atoms[i].label)
        neighbours = list(species.graph.neighbors(i))

        if len(neighbours) <= max_valance:
            # All is well
            continue

        # Get the atom indexes sorted by the closest to atom i
        closest_atoms = sorted(
            neighbours, key=lambda j: species.get_distance(i, j))

        # Delete all the bonds to atom(s) j that are above the maximal valance
        for j in closest_atoms[max_valance:]:
            species.graph.remove_edge(i, j)
    return None
