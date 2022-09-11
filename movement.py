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
        self.moved = False
        self.move()

    def move(self):
        """
        Compute movement.
        """
        if self.type == "end":
            self.end_movement()
        elif self.type == "corner":
            self.corner_movement()
        elif self.type == "crankshaft":
            self.crankshaft_movement()
        elif self.type == "pull":
            self.pull_movement()

    def end_movement(self):
        """
        Compute end movement.

        Returns
        -------
        tuple of coordinates of the new position of the residue
        """
        neighbor_residue = self.lattice.protein.get_consecutive(self.residue)

        # end residues have exactly one neighbor
        if len(neighbor_residue) == 1:
            empty_neighbors = self.lattice.empty_neighbors(
                neighbor_residue[0].get_coords())

            # if another position is available
            if empty_neighbors:
                random_neighbor = empty_neighbors[np.random.choice(
                    len(empty_neighbors))]
                self.lattice.move_residue(self.residue, random_neighbor)
                self.moved = True

    def corner_movement(self):
        """
        Compute corner movement.

        Returns
        -------
        tuple of coordinates of the new position of the residue
        """
        neighbors_residues = self.lattice.protein.get_consecutive(self.residue)

        # corner residues have exactly two neighbors
        if len(neighbors_residues) == 2:
            neighbors_residues_coords = tuple(
                res.get_coords() for res in neighbors_residues
            )

            corner_position = tuple(
                abs(i + j) - k
                for i, j, k in zip(*neighbors_residues_coords, self.residue.get_coords())
            )

            # if the corner position is available
            if self.lattice.is_empty(corner_position):
                self.lattice.move_residue(self.residue, corner_position)
                self.moved = True

    def crankshaft_movement(self):
        """
        Compute crankshaft movement.

        Returns
        -------
        tuple of coordinates of the new position of the residues
        """
        start_index = self.residue.index
        neighbors_residues = self.lattice.protein.get_consecutive(self.residue)

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
                                    self.residue.get_coords(),
                                )
                            )

                            # other residue to move
                            other_residue = self.lattice.protein.get_residue(
                                position[2]
                            )
                            # new position of the other residue
                            new_position_j = tuple(
                                (2 * (i - j)) + j
                                for i, j in zip(
                                    corner_candidates_coords[1],
                                    other_residue.get_coords(),
                                )
                            )

                            # if the new positions are available
                            if self.lattice.is_empty(new_position_i) \
                                    and self.lattice.is_empty(new_position_j):
                                self.lattice.move_residue(self.residue,
                                                          new_position_i)
                                self.lattice.move_residue(other_residue,
                                                          new_position_j)
                                moved = True
                                return

    def pull_movement(self):
        """
        Compute pull movement.

        Returns
        -------
        tuple of coordinates of the new position of the residues
        """
        neighbors_residues = self.lattice.protein.get_consecutive(self.residue)

        # pull residues have exactly two neighbors
        if len(neighbors_residues) == 2:
            for index, neighbor_plus_1 in enumerate(neighbors_residues):
                # check the direction between start residue and neighbor
                # invert that direction to get the new position of the start residue
                side_direction = tuple(
                    map(lambda a, b: 1 - abs(a - b),
                        self.residue.get_coords(),
                        neighbor_plus_1.get_coords())
                )

                # side positions to move to
                c = tuple(map(add, self.residue.get_coords(), side_direction))
                l = tuple(map(add, neighbor_plus_1.get_coords(), side_direction))

                if self.lattice.is_empty(c) and self.lattice.is_empty(l):
                    neighbor_minus_1 = neighbors_residues[1 - index]

                    # save indexes of the residues to move
                    i_minus_2 = self.residue.get_coords()
                    i_minus_1 = neighbor_minus_1.get_coords()

                    self.lattice.move_residue(self.residue, l)
                    self.lattice.move_residue(neighbor_minus_1, c)

                    # which way to head in the protein, either 1 or -1
                    direction = self.residue.index - neighbor_minus_1.index
                    next_residue = self.lattice.protein.get_residue(
                        neighbor_minus_1.index + direction)

                    # loop over residues in that direction
                    while next_residue is not None:
                        # redefine the new positions
                        new_coords = i_minus_2
                        i_minus_2 = i_minus_1
                        i_minus_1 = next_residue.get_coords()

                        # move the residue
                        self.lattice.move_residue(next_residue, new_coords)

                        # check if the conformation is valid
                        if self.lattice.is_valid():
                            moved = True
                            return
                        else:
                            next_residue = self.lattice.protein.get_residue(
                                next_residue.index + direction)

    def __str__(self):
        return f"{self.type} movement of {self.residue}"
