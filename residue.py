"""Create and manipulate residues.

Simple class for handling HP residues. HP residues are a simplification of
protein residues, classified as either hydrophobic (H) or polar (P).
"""


class Residue(object):
    """
    Create and manipulate residues.

    Attributes
    ----------
    index : int
        index of the residu in the sequence
    typeHP : char
        type of the residu ; Hydrophobic (H) or Polar (P)
    coordX : int
        x coordinate of the residue
    coordY : int
        y coordinate of the residue

    Methods
    -------
    get_coords():
        Return the coordinates of the residue.
    set_coords(tuple):
        Set the coordinates of the residue.
    """

    def __init__(self, input_index, input_type):
        """
        create Residue object

        Parameters
        ----------
        input_index : int
            index of the residue in the input sequence
        input_type : char
            residue type (H or P)
        """
        self.index = input_index
        self.typeHP = input_type
        self.coordX = None
        self.coordY = None

    def get_coords(self):
        """
        Get the coordinates of the residue

        Returns
        -------
        tuple of residue's coordinates
        """
        return(self.coordX, self.coordY)

    def set_coords(self, input_coords):
        """
        Set the coordinates of the residue

        Parameters
        ----------
        input_coords : tuple
            tuple of coordinates
        """
        self.coordX, self.coordY = input_coords

    def is_consecutive(self, residue2):
        """Check if two residues are consecutive.

        Parameters
        ----------
        residue2 : Residue
            Residue to check.

        Returns
        -------
        True if the residues are consecutive, False otherwise.
        """
        return abs(self.index - residue2.index) == 1

    def __str__(self):
        return self.typeHP