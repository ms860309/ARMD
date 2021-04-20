# third party
import numpy as np


class Species:
    def __init__(self, atoms):
        self.atoms = atoms
        self.n_atoms = 0 if atoms is None else len(atoms)
        self.graph = None

    def get_coordinates(self):
        """Return a np.ndarray of size n_atoms x 3 containing the xyz
        coordinates of the molecule"""
        return np.array([atom.coord for atom in self.atoms])
