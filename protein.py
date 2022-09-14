"""Create and manipulate proteins.

Simple class for handling proteins as chains of HP residues.
See the `Residue` class for more information on HP residues.
"""

# standard library
from copy import deepcopy
import math
from sqlite3 import connect

# local
from residue import Residue


class Protein():
    """
    Create and manipulate proteins as chains of residues.

    Attributes
    ----------
    length : int
        Length of the protein.
    sequence : str
        Sequence of residues.
    residues : list
        List of residues.

    Methods
    -------
    get_residue(index):
        Return the residue at the given index.
    get_H_residues():
        Return the list of hydrophobic residues.
    get_consecutive(residue):
        Return the neighbors of the given residue.
    is_end(residue):
        Check if the residue is at either end of the protein.
    is_corner(residue):
        Check if the residue is in a corner of the protein.
    """

    def __init__(self, sequence):
        """
        create Protein object

        Parameters
        ----------
        sequence : str
            sequence of residues
        """
        self.length = len(sequence)
        self.sequence = sequence
        self.residues = [Residue(i, sequence[i]) for i in range(self.length)]

    def get_residue(self, index):
        """
        Get the residue at the given index.

        Parameters
        ----------
        index : int
            index of the residue in the sequence

        Returns
        -------
        residue at the given index
        """
        if index >= self.length or index < 0:
            return None
        return self.residues[index]

    def get_H_residues(self):
        """
        Return the list of hydrophobic residues.

        Returns
        -------
        List of hydrophobic residues.
        """
        return [residue for residue in self.residues if residue.typeHP == "H"]

    def get_consecutive(self, residue):
        """
        Return the neighbors of the given residue.

        Parameters
        ----------
        residue : Residue
            Residue to check.

        Returns
        -------
        List of neighbors.
        """
        if residue.index == 0:
            return [self.residues[1]]
        elif residue.index == self.length - 1:
            return [self.residues[-2]]
        else:
            return [self.residues[residue.index - 1], self.residues[residue.index + 1]]

    def is_end(self, residue):
        """Check if the residue is at the end of the protein.

        Parameters
        ----------
        residue : Residue
            Residue to check.

        Returns
        -------
        True if the residue is at the end of the protein, False otherwise.
        """
        return residue.index in (0, self.length - 1)

    def is_corner(self, residue):
        """Check if the residue is in a corner of the protein.

        Parameters
        ----------
        residue : Residue
            Residue to check.

        Returns
        -------
        True if the residue is in a corner of the protein, False otherwise.
        """
        neighbors_residues = self.get_consecutive(residue)

        # corner residues have exactly two neighbors
        if len(neighbors_residues) == 2:
            neighbors_residues_coords = [
                res.get_coords() for res in neighbors_residues]

            # check that the two neighbors form a corner
            if math.prod(map(lambda a, b: abs(a - b), *neighbors_residues_coords)) == 1:
                return True
        return False

    def write_pdb(self, filename):
        """
        Write a PDB file for the protein.

        Parameters
        ----------
        filename : str
            Name of the PDB file to write.
        """
        with open(filename, "w") as handle:
            atom = ''
            connect = ''
            for i, residue in enumerate(self.residues):
                atom += (f"ATOM  {i:>5}{'C' if residue.typeHP == 'H' else 'O':>4}"
                         f" ACY A {i:>3} {residue.coordI:8.3f}{residue.coordJ:8.3f}{0:8.3f}"
                         f"  1.00  0.00           {'C' if residue.typeHP == 'H' else 'O'}  \n"
                         )
                connect += (f"CONECT{i:>5}{((i+1) % self.length):>5}\n")
            handle.write(atom)
            handle.write(connect)
            handle.write("END\n")

    def __str__(self):
        return self.sequence

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for key, value in self.__dict__.items():
            setattr(result, key, deepcopy(value, memo))
        return result
