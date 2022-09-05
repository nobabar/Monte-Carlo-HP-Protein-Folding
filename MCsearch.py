import copy
import numpy as np

import vshd


def MCsearch(nsteps, lattice_input, temperature):
    """
    Perform a Monte Carlo search of the lattice.

    Parameters
    ----------
    nsteps : int
        Number of steps to perform.
    lattice : Lattice
        Lattice on which to perform the search.
    temperature : float
        Temperature of the search.

    Returns
    -------
    Lattice
        Lattice with the protein placed on it.
    """
    # copy the lattice
    lattice = copy.deepcopy(lattice_input)

    # compute the initial energy
    energy = lattice.calculate_energy()
    print(f"Initial energy: {energy}")
    new_energy = energy

    # perform the search
    for i in range(nsteps):
        # choose a random residue
        residue = np.random.choice(lattice.protein.residues)

        # compute the movement
        movement = None
        if lattice.protein.is_end(residue):
            movement = vshd.end_movement(lattice, residue)
        else:
            movement = vshd.corner_movement(lattice, residue)

        # if a movement is possible
        if movement:
            if residue.typeHP == "H":
                # compute the new energy
                new_energy = energy + \
                    lattice.calculate_energy_change(residue, movement)

            # if the new energy is lower or if the Boltzmann condition is met
            if new_energy <= energy or np.random.random() < np.exp(-(new_energy - energy) / temperature):
                # move the residue
                lattice.move_residue(residue, movement)
                print(f"Iteration {i}, energy : {new_energy}")
                lattice.draw_grid()

                # update the energy
                energy = new_energy

    return lattice
