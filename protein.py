"""Create and manipulate proteins.

Simple class for handling proteins as chains of HP residues.
See the `Residue` class for more information on HP residues.
"""

from residue import Residue


class Protein(object):
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
            raise IndexError("Index out of range")
        return self.residues[index]

    def get_H_residues(self):
        """
        Return the list of hydrophobic residues.

        Returns
        -------
        List of hydrophobic residues.
        """
        return [residue for residue in self.residues if residue.typeHP == "H"]

    def get_neighbors(self, residue):
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
            return [self.residues[residue.index - 1],
                    self.residues[residue.index + 1]]

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
        return residue.index == 0 or residue.index == self.length - 1

    def __str__(self):
        return self.sequence
