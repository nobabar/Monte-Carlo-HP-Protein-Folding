import copy
import numpy as np
from MCsearch import MCsearch

from movement import Movement


def REMCsearch(
    n_replica, energy_cutoff, max_steps, local_steps, temperature_min, temperature_max, lattice_input
):
    """
    Perform a Replica Exchange Monte Carlo search on the lattice.

    Parameters
    ----------
    n_replica : int
        Number of replicas to use.
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
    lattices = [copy.deepcopy(lattice_input) for _ in range(n_replica)]
    temperatures = np.linspace(temperature_min, temperature_max, n_replica)
    offset = 0
    energy = 0
    step = 0
    while energy > energy_cutoff and step < max_steps:
        for replica in range(n_replica):
            lattice = MCsearch(
                local_steps, temperatures[replica], lattices[replica])
            if lattice.calculate_energy() < lattices[replica].calculate_energy():
                lattices[replica] = lattice

        # if the replica with the minimum energy pass the cutoff
        if min([l.calculate_energy() for l in lattices]) < energy_cutoff:
            break

        i = offset
        while i < (n_replica - 1):
            j = i + 1

            # Boltzmann constant
            K_b = 0.0019872041

            # product of the energy difference and inverse temperature difference
            delta = ((1 / (temperatures[j] * K_b)) - (1 / (temperatures[i] * K_b))) * (
                lattices[i].calculate_energy() - lattices[j].calculate_energy())

            if delta <= 0 or np.random.random() <= np.exp(-delta):
                lattices[i], lattices[j] = lattices[j], lattices[i]
            i += 2
        offset = 1 - offset
        step += 1
    # return the lattice with the lowest energy
    return min(lattices, key=lambda x: x.calculate_energy())
