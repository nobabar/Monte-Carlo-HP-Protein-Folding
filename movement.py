import copy
import numpy as np
from operator import add


class Movement(object):
    """
    Move residue in the lattice.

    Attributes
    ----------
    type : str
        Type of movement.
    lattice : Lattice
        Lattice in which to make the movement.
    residue : Residue
        Residue to move.
    destination : tuple
        destination of the residue.

    Methods
    -------
    move()
        Move residue in the lattice.
    end_movement()
        Compute end movement.
    corner_movement()
        Compute corner movement.
    crankshaft_movement()
        Compute crankshaft movement.
    pull_movement()
        Compute pull movement.
    """

    def __init__(self, type, lattice, residue):
        """
        Initialize a movement.

        Parameters
        ----------
        type : str
            Type of movement.
        lattice : Lattice
            Lattice in which to make the movement.
        residue : Residue
            Residue to move.
        """
        self.type = type
        self.lattice = copy.deepcopy(lattice)
        self.residue = residue
        self.destination = None
        self.move()

    def move(self):
        """
        Compute movement.
        """
        if self.type == "end":
            self.destination = self.end_movement()
        elif self.type == "corner":
            self.destination = self.corner_movement()
        elif self.type == "crankshaft":
            self.destination = self.crankshaft_movement()
        elif self.type == "pull":
            self.destination = self.pull_movement()

    def end_movement(self):
        """
        Compute end movement.

        Returns
        -------
        tuple of coordinates of the new position of the residue
        """
        neighbor_residue = self.lattice.protein.get_consecutive(self.residue)
        empty_neighbors = self.lattice.empty_neighbors(neighbor_residue.get_coords())

        # if another position is available
        if empty_neighbors:
            random_neighbor = np.random.choice(len(empty_neighbors))

            if self.residue.typeHP == "H":
                # compute the new energy
                new_energy = new_energy + self.lattice.calculate_energy_change(
                    self.residue, self.destination
                )

            return empty_neighbors[random_neighbor]
        return

    def corner_movement(self):
        """
        Compute corner movement.

        Returns
        -------
        tuple of coordinates of the new position of the residue
        """
        neighbors_residues = self.lattice.protein.get_consecutive(self.residue[0])
        neighbors_residues_coords = tuple(
            res.get_coords() for res in neighbors_residues
        )

        corner_position = tuple(
            abs(i + j) - k
            for i, j, k in zip(*neighbors_residues_coords, self.residue[0].get_coords())
        )

        # if the corner position is available
        if self.lattice.is_empty(corner_position):
            return corner_position
        return

    def crankshaft_movement(self):
        """
        Compute crankshaft movement.

        Returns
        -------
        tuple of coordinates of the new position of the residues
        """
        start_index = self.residue[0].index
        neighbors_residues = self.lattice.protein.get_consecutive(self.residue[0])

        # crankshaft residues have exactly two neighbors
        if len(neighbors_residues) == 2:
            # check that i-1 and i+2 are corner residues
            # or that i-2 and i+1 are corner residues
            for position in [
                (start_index - 1, start_index + 2, start_index + 1),
                (start_index + 1, start_index - 2, start_index - 1),
            ]:
                corner_candidates = [
                    self.lattice.protein.get_residue(position[0]),
                    self.lattice.protein.get_residue(position[1]),
                ]

                if not None in corner_candidates:
                    if all(
                        [
                            self.lattice.protein.is_corner(candidate)
                            for candidate in corner_candidates
                        ]
                    ):
                        corner_candidates_coords = tuple(
                            res.get_coords() for res in corner_candidates
                        )

                        # check that they also are neighbors
                        if self.lattice.are_neighbors(*corner_candidates_coords):
                            # new position of initial residue
                            new_position_i = tuple(
                                (2 * (i - j)) + j
                                for i, j in zip(
                                    corner_candidates_coords[0],
                                    self.residue[0].get_coords(),
                                )
                            )

                            # other residue to move
                            other_residue = self.lattice.protein.get_residue(
                                position[2]
                            )
                            self.residue = [self.residue[0], other_residue]
                            # new position of the other residue
                            new_position_j = tuple(
                                (2 * (i - j)) + j
                                for i, j in zip(
                                    corner_candidates_coords[1],
                                    other_residue.get_coords(),
                                )
                            )

                            # if the new positions are available
                            if self.lattice.is_empty(
                                new_position_i
                            ) and self.lattice.is_empty(new_position_j):
                                return [new_position_i, new_position_j]
        return

    def pull_movement(self):
        """
        Compute pull movement.

        Returns
        -------
        tuple of coordinates of the new position of the residues
        """
        start_index = self.residue[0].index
        neighbors_residues = self.lattice.protein.get_consecutive(self.residue[0])

        # pull residues have exactly two neighbors
        if len(neighbors_residues) == 2:
            for neighbor in neighbors_residues:
                # check the direction between start residue and neighbor
                # invert that direction to get the new position of the start residue
                side_direction = tuple(
                    map(lambda a, b: 1 - abs(a - b), self.residue[0], neighbor)
                )

                # side positions to move to
                c = tuple(map(add, self.residue[0], side_direction))
                l = tuple(map(add, neighbor, side_direction))

                if self.lattice.is_empty(c) and self.lattice.is_empty(l):
                    ...

        return

    def __str__(self):
        return f"{self.type} movement of {self.residue} from {[res.get_coords() for res in self.residue]} to {self.destination}"
