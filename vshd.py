import numpy as np


class VSHD(object):
    """
    Move residue in the lattice.

    Attributes
    ----------
    lattice : Lattice
        Lattice in which to make the movement.
    residues : Residue
        Residue(s) to move.
    destinations : tuple
        destination(s) of the residue(s).

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
    """
    def __init__(self, type, lattice, residue):
        """
        Initialize a VHSD movement.
        
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
        self.lattice = lattice
        self.residues = [residue]
        self.destinations = []
        self.move()

    def move(self):
        """
        Compute movement.

        Parameters
        ----------
        lattice : Lattice
            Lattice in which to make the movement.
        residue : Residue
            Residue to move.
        """
        if self.type == "end":
            self.destinations.append(self.end_movement())
        elif self.type == "corner":
            self.destinations.append(self.corner_movement())
        elif self.type == "crankshaft":
            self.destinations.append(self.crankshaft_movement())

    def end_movement(self):
        """
        Compute end movement.

        Parameters
        ----------
        lattice : Lattice
            Lattice in which to make the movement.
        residue : Residue
            End residue to move.
        """
        neighbor_residue = self.lattice.protein.get_consecutive(self.residues[0])
        empty_neighbors = self.lattice.empty_neighbors(neighbor_residue[0].get_coords())

        # if another position is available
        if empty_neighbors:
            random_neighbor = np.random.choice(len(empty_neighbors))
            return empty_neighbors[random_neighbor]
        return


    def corner_movement(self):
        """
        Compute corner movement.

        Parameters
        ----------
        lattice : Lattice
            Lattice in which to make the movement.
        residue : Residue
            Corner residue to move.
        """
        neighbors_residues = self.lattice.protein.get_consecutive(self.residues[0])
        neighbors_residues_coords = tuple(
            res.get_coords() for res in neighbors_residues)

        corner_position = tuple(abs(
            i + j) - k for i, j, k in zip(*neighbors_residues_coords,
                                        self.residues[0].get_coords()))

        # if the corner position is available
        if self.lattice.is_empty(corner_position):
            return corner_position
        return


    def crankshaft_movement(self):
        """
        Compute crankshaft movement.

        Parameters
        ----------
        lattice : Lattice
            Lattice in which to make the movement.
        residue : Residue
            Residue to move.
        """
        start_index = self.residues[0].index
        neighbors_residues = self.lattice.protein.get_consecutive(self.residues[0])

        # crankshaft residues have exactly two neighbors
        if len(neighbors_residues) == 2:
            # check that i-1 and i+2 are corner residues
            # or that i-2 and i+1 are corner residues
            for position in [(start_index - 1, start_index + 2, start_index + 1),
                            (start_index + 1, start_index - 2, start_index - 1)]:
                corner_candidates = [self.lattice.protein.get_residue(position[0]),
                                    self.lattice.protein.get_residue(position[1])]

                if not None in corner_candidates:
                    if all([self.lattice.protein.is_corner(candidate) for candidate in corner_candidates]):
                        corner_candidates_coords = tuple(
                            res.get_coords() for res in corner_candidates)

                        # check that they also are neighbors
                        if self.lattice.are_neighbors(*corner_candidates_coords):
                            # new position of initial residue
                            new_position_i = tuple(
                                (2 * (i - j)) + j for i, j in zip(corner_candidates_coords[0],
                                                                self.residues[0].get_coords()))

                            # other residue to move
                            other_residue = self.lattice.protein.get_residue(position[2])
                            self.residues = [self.residues[0], other_residue]
                            # new position of the other residue
                            new_position_j = tuple(
                                (2 * (i - j)) + j for i, j in zip(corner_candidates_coords[1],
                                                                other_residue.get_coords()))

                            # if the new positions are available
                            if self.lattice.is_empty(new_position_i) and self.lattice.is_empty(new_position_j):
                                print("crankshaft!")
                                return [new_position_i, new_position_j]
        return
    
    def __str__(self):
        return f"{self.type} movement of {self.residues} from {[res.get_coords() for res in self.residues]} to {self.destinations}"
