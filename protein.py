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
        return self.residues[index]

    def __str__(self):
        return self.sequence
