import numpy as np
import random

from residue import Residue
from protein import Protein


class Lattice(object):
    """
    Lattice grid to restrict residue placement.

    Attributes
    ----------
    size : int
        Size of the grid.
    grid : numpy.ndarray
        Grid of characters.
    protein : Protein
        Protein to place on the grid.

    Methods
    -------
    place_residue(residue, coords):
        Place a residue on the grid.
    remove_residue(coords):
        Remove a residue from the grid.
    fill_grid(protein):
        Fill the grid with a random generated conformation.
    is_empty(coords):
        Check if the given coordinates are empty.
    empty_neighbors(coords):
        Return the empty neighbors of the given coordinates.
    draw_grid():
        Draw the lattice in the terminal.
    """

    def __init__(self, protein):
        """
        Initialize an empty lattice.
        """
        self.size = protein.length
        self.grid = np.ndarray(shape=(self.size, self.size), dtype=Residue)
        self.protein = protein
        self.fill_grid()

    def place_residue(self, residue, coords):
        """
        Place a residue on the grid.

        Parameters
        ----------
        residue : Residue
            Residue to place on the grid.
        coords : tuple
            Coordinates of the residue.
        """
        self.grid[coords] = residue
        residue.set_coords(coords)

    def remove_residue(self, coords):
        """
        Remove a residue from the grid.

        Parameters
        ----------
        coords : tuple
            Coordinates of the residue to remove.
        """
        self.grid[coords] = None

    def fill_grid(self):
        """
        Fill the grid with a random generated conformation.

        Parameters
        ----------
        protein : Protein
            Protein to place on the grid.
        """
        # place the first residue in the middle of the grid
        midpoint = (int(self.size // 2), int(self.size // 2))
        self.place_residue(self.protein.get_residue(0), midpoint)

        # place the remaining residues
        coords = midpoint
        dead_ends = []
        for i in range(1, self.protein.length):
            res = self.protein.get_residue(i)
            empty_neighbors = self.empty_neighbors(coords)

            for dead_end in dead_ends:
                empty_neighbors.remove(
                    dead_end) if dead_end in empty_neighbors else None

            if len(empty_neighbors) == 0:
                dead_ends.append(coords)
                self.remove_residue(coords)
                coords = self.protein.get_residue(i - 1).get_coords()
                i -= 1
            else:
                random_neighbor = random.choice(empty_neighbors)
                self.place_residue(res, random_neighbor)
                coords = random_neighbor

    def is_empty(self, coords):
        """
        Check if the given coordinates are empty.

        Parameters
        ----------
        coords : tuple
            Coordinates to check.

        Returns
        -------
        True if the coordinates are empty, False otherwise.
        """
        return self.grid[coords] is None

    def empty_neighbors(self, coords):
        """
        Return the empty neighbors of the given coordinates.

        Parameters
        ----------
        coords : tuple
            Coordinates to check.

        Returns
        -------
        List of empty neighbors.
        """
        neighbors = []
        if coords[0] > 0:
            neighbors.append((coords[0] - 1, coords[1]))
        if coords[0] < self.size - 1:
            neighbors.append((coords[0] + 1, coords[1]))
        if coords[1] > 0:
            neighbors.append((coords[0], coords[1] - 1))
        if coords[1] < self.size - 1:
            neighbors.append((coords[0], coords[1] + 1))
        return [n for n in neighbors if self.is_empty(n)]

    def draw_grid(self):
        """
        Draw the lattice in the terminal
        """
        for i in range(self.size):
            print("+---" * self.size + "+")
            for j in range(self.size):
                if self.grid[i, j] is None:
                    print("|   ", end="")
                else:
                    print("|" + f"{str(self.grid[i, j]):^3}", end="")
            print("|")
        print("+---" * self.size + "+")
