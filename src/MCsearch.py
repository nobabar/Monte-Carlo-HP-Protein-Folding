import copy
import numpy as np

from src.movement import Movement


def MCsearch(n_steps, temperature, lattice_input):
    """
    Perform a Monte Carlo search of the lattice.

    Parameters
    ----------
    n_steps : int
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
    new_energy = energy

    # perform the search
    for _ in range(n_steps):
        # choose a random residue
        residue = np.random.choice(lattice.protein.residues)

        # compute the movement
        movements = []
        if lattice.protein.is_end(residue):
            movements.append(Movement("end", lattice, residue))
        else:
            if lattice.protein.is_corner(residue):
                movements.append(Movement("corner", lattice, residue))
                movements.append(Movement("crankshaft", lattice, residue))
            movements.append(Movement("pull", lattice, residue))

        # filter movements
        movements = [m for m in movements if m.moved]

        if movements:
            random_movement = np.random.choice(movements)
            new_energy = random_movement.lattice.calculate_energy()

            # Boltzmann constant
            K_b = 0.0019872041

            # if the new energy is lower or if the Boltzmann condition is met
            if new_energy <= energy or np.random.random() < np.exp(
                -(new_energy - energy) / (temperature * K_b)
            ):
                # update the lattice
                lattice = random_movement.lattice

                # print(f"Iteration {i}, energy : {new_energy}")
                # lattice.draw_grid()

                # update the energy
                energy = new_energy
    return lattice
