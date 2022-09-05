import numpy as np

from residue import Residue


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

    def __init__(self, protein, initial_placement_mode="linear"):
        """
        Initialize an empty lattice.

        Parameters
        ----------
        protein : Protein
            Protein to place on the grid.
        initial_placement_mode : str
            Initial placement mode of the protein on the grid.
            Either linear or random.
        """
        self.size = protein.length * 2
        self.grid = np.ndarray(shape=(self.size, self.size), dtype=Residue)
        self.protein = protein
        self.fill_grid(initial_placement_mode)

    def get_residue(self, coords):
        """
        Return the residue at the given coordinates.

        Parameters
        ----------
        coords : tuple
            Coordinates of the residue.

        Returns
        -------
        Residue at the given coordinates.
        """
        return self.grid[coords]

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

    def move_residue(self, residue, coords):
        """
        Move a residue to the given coordinates.

        Parameters
        ----------
        residue : Residue
            Residue to move.
        coords : tuple
            Coordinates to move the residue to.
        """
        self.remove_residue(residue.get_coords())
        self.place_residue(residue, coords)

    def remove_residue(self, coords):
        """
        Remove a residue from the grid.

        Parameters
        ----------
        coords : tuple
            Coordinates of the residue to remove.
        """
        self.grid[coords] = None

    def fill_grid(self, mode):
        """
        Fill the grid with a random generated conformation.

        Parameters
        ----------
        protein : Protein
            Protein to place on the grid.
        """
        if mode == "linear":
            i = j = self.size // 2 - self.protein.length // 2
            start_point = (i, j)
            self.place_residue(self.protein.get_residue(0), start_point)

            for index in range(1, self.protein.length):
                self.place_residue(
                    self.protein.get_residue(index), (i, j + index))

        elif mode == "random":
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
                    random_neighbor = empty_neighbors[np.random.choice(
                        len(empty_neighbors))]
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
        Return the empty neighbors around a given coordinates.

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

    def occupied_neighbors(self, coords):
        """
        Return the occupied neighbors around a given coordinates.

        Parameters
        ----------
        coords : tuple
            Coordinates to check.

        Returns
        -------
        List of occupied neighbors.
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
        return [self.get_residue(n) for n in neighbors if not self.is_empty(n)]

    def calculate_energy(self):
        """
        Calculate the energy of the lattice.

        Returns
        -------
        Energy of the lattice.
        """
        energy = 0

        h_residues = self.protein.get_H_residues()
        for res in h_residues:
            neighbors = self.occupied_neighbors(res.get_coords())
            for neighbor in neighbors:
                if not res.is_consecutive(neighbor) and neighbor.typeHP == "H":
                    energy -= 1
        # don't count each bound twice
        return energy / 2

    def calculate_energy_change(self, residue, movement):
        """
        Calculate the energy change of a residue after a movement.

        Parameters
        ----------
        residue : Residue
            Residue to check.
        movement : tuple
            Movement to apply.

        Returns
        -------
        Energy change of the residue after the movement.
        """
        energy_change = 0

        initial_neighbors = self.occupied_neighbors(residue.get_coords())
        for neighbor in initial_neighbors:
            if not residue.is_consecutive(neighbor) and neighbor.typeHP == "H":
                energy_change += 1

        new_neighbors = self.occupied_neighbors(movement)
        for neighbor in new_neighbors:
            if not residue.is_consecutive(neighbor) and neighbor.typeHP == "H":
                energy_change -= 1
        return energy_change

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
