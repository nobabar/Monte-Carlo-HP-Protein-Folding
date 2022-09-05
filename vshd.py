import numpy as np
import copy
import math


class VHSD(object):
    """
    Class to compute VHSD mouvements

    Attributes
    ----------
    lattice : Lattice
        Lattice containing the protein.
    protein : Protein
        Protein to move.
    """

    def __init__(self, input_lattice, input_protein):
        """
        Initialize the VHSD mouvements

        Parameters
        ----------
        input_lattice : Lattice
            Lattice containing the protein.
        input_protein : Protein
            Protein to move.
        """
        self.lattice = copy.deepcopy(input_lattice)
        self.protein = copy.deepcopy(input_protein)

    def end_movement(self, input_residue):
        """
        Compute end movement.

        Parameters
        ----------
        input_residue : Residue
            End residue to move.
        """
        # don't manipulate the input_residu but own copy
        residue = self.protein.get_residue(input_residue.index)

        neighbor_residue = self.protein.get_neighbors(residue)
        if neighbor_residue.length != 1:
            raise ValueError(
                "Residue isn't an end residue as it has not exactly one neighbor")

        empty_neighbors = self.grid.empty_neighbors(neighbor_residue)
        random_neighbor = np.random.choice(empty_neighbors)
        self.lattice.place_residue(residue, random_neighbor)

    def corner_movement(self, input_residue):
        """
        Compute corner movement.

        Parameters
        ----------
        input_residue : Residue
            Corner residue to move.
        """
        # don't manipulate the input_residu but own copy
        residue = self.protein.get_residue(input_residue.index)

        neighbor_residue = self.protein.get_neighbors(residue)
        if neighbor_residue.length != 2:
            raise ValueError(
                "Residue isn't a corner residue as it has not exactly two neighbors")

        math.prod([abs(i - j) for i, j in zip(c, b)])
