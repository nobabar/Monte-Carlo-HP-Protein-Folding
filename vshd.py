import numpy as np
import math


def end_movement(lattice, residue):
    """
    Compute end movement.

    Parameters
    ----------
    lattice : Lattice
        Lattice in which to make the movement.
    residue : Residue
        End residue to move.
    """
    neighbor_residue = lattice.protein.get_neighbors(residue)
    if len(neighbor_residue) != 1:
        raise ValueError(
            "Residue isn't an end residue as it has not exactly one neighbor")

    empty_neighbors = lattice.empty_neighbors(neighbor_residue[0].get_coords())

    # if another position is available
    if empty_neighbors:
        random_neighbor = np.random.choice(len(empty_neighbors))
        return empty_neighbors[random_neighbor]
    else:
        return None


def corner_movement(lattice, residue):
    """
    Compute corner movement.

    Parameters
    ----------
    lattice : Lattice
        Lattice in which to make the movement.
    residue : Residue
        Corner residue to move.
    """
    neighbors_residues = lattice.protein.get_neighbors(residue)

    # corner residues have exactly two neighbors
    if len(neighbors_residues) == 2:
        neighbors_residues_coords = tuple(
            res.get_coords() for res in neighbors_residues)

        # check that the two neighbors form a corner
        if math.prod([abs(i - j) for i, j in zip(*neighbors_residues_coords)]):
            corner_position = tuple(abs(
                i + j) - k for i, j, k in zip(*neighbors_residues_coords, residue.get_coords()))

            # if the corner position is available
            if lattice.is_empty(corner_position):
                return corner_position
    return None
