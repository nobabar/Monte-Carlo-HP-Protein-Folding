import numpy as np


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
    neighbor_residue = lattice.protein.get_consecutive(residue)
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
    neighbors_residues = lattice.protein.get_consecutive(residue)
    neighbors_residues_coords = tuple(
        res.get_coords() for res in neighbors_residues)

    corner_position = tuple(abs(
        i + j) - k for i, j, k in zip(*neighbors_residues_coords,
                                      residue.get_coords()))

    # if the corner position is available
    if lattice.is_empty(corner_position):
        return corner_position
    return None


def crankshaft_movement(lattice, residue):
    """
    Compute crankshaft movement.

    Parameters
    ----------
    lattice : Lattice
        Lattice in which to make the movement.
    residue : Residue
        Residue to move.
    """
    start_index = residue.index
    neighbors_residues = lattice.protein.get_consecutive(residue)

    # crankshaft residues have exactly two neighbors
    if len(neighbors_residues) == 2:
        neighbors_residues_coords = tuple(
            res.get_coords() for res in neighbors_residues)

        # check that i-1 and i+2 are corner residues
        # or that i-2 and i+1 are corner residues
        for position in [(start_index - 1, start_index + 2, start_index + 1),
                         (start_index + 1, start_index - 2, start_index - 1)]:
            corner_candidates = [lattice.protein.get_residue(position[0]),
                                 lattice.protein.get_residue(position[1])]

            if all([lattice.protein.is_corner(candidate) for candidate in corner_candidates]):
                corner_candidates_coords = tuple(
                    res.get_coords() for res in corner_candidates)

                # check that they also are neighbors
                if lattice.are_neighbors(*corner_candidates_coords):
                    # new position of initial residue
                    new_position_i = tuple(
                        (2 * (i - j)) + j for i, j in zip(corner_candidates_coords[0],
                                                          residue.get_coords()))

                    # other residue to move
                    other_residue = lattice.protein.get_residue(position[2])
                    # new position of the other residue
                    new_position_j = tuple(
                        (2 * (i - j)) + j for i, j in zip(corner_candidates_coords[1],
                                                          other_residue.get_coords()))

                    # if the new positions are available
                    if lattice.is_empty(new_position_i) and lattice.is_empty(new_position_j):
                        return new_position_i, new_position_j
    return None
