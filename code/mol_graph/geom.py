# third party
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import distance_matrix


def get_neighbour_list(species, atom_i):
    """Calculate a neighbour list from atom i as a list of atom labels
    Arguments:
        atom_i (int): index of the atom
        species (autode.species.Species):
    Returns:
        (list(int)): list of atom ids in ascending distance away from atom_i
    """
    coords = species.get_coordinates()
    distance_vector = cdist(np.array([coords[atom_i]]), coords)[0]

    dists_and_atom_labels = {}
    for atom_j, dist in enumerate(distance_vector):
        dists_and_atom_labels[dist] = species.atoms[atom_j].label

    atom_label_neighbour_list = []
    for dist, atom_label in sorted(dists_and_atom_labels.items()):
        atom_label_neighbour_list.append(atom_label)

    return atom_label_neighbour_list


def get_points_on_sphere(n_points, r=1):
    """
    Find n evenly spaced points on a sphere using the "How to generate
    equidistributed points on the surface of a sphere" by Markus Deserno, 2004.
    Arguments:
        n_points (int): number of points to generate
        r (float): radius of the sphere
    Returns:
        (list(np.ndarray))
    """
    points = []

    a = 4.0 * np.pi * r**2 / n_points
    d = np.sqrt(a)
    m_theta = int(np.round(np.pi / d))
    d_theta = np.pi / m_theta
    d_phi = a / d_theta

    for m in range(m_theta):
        theta = np.pi * (m + 0.5) / m_theta
        m_phi = int(np.round(2.0 * np.pi * np.sin(theta)/d_phi))

        for n in range(m_phi):
            phi = 2.0 * np.pi * n / m_phi
            point = np.array([r * np.sin(theta) * np.cos(phi),
                              r * np.sin(theta) * np.sin(phi),
                              r * np.cos(theta)])

            points.append(point)

    return points
