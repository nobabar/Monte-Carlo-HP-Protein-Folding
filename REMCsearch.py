import copy
import numpy as np
from MCsearch import MCsearch

from movement import Movement


def REMCsearch(n_replica, n_steps, n_local_steps, t_min, t_max, lattice_input):
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
    temperatures = np.linspace(t_min, t_max, n_replica)
    offset = 0
    for _ in range(n_steps):
        for replica in range(n_replica):
            lattice = MCsearch(n_local_steps, lattices[replica], temperatures[replica])
            if lattice.calculate_energy() < lattices[replica].calculate_energy():
                lattices[replica] = lattice

        i = offset + 1
        while i < (n_replica - 1):
            j = i + 1

            # product of the energy difference and inverse temperature difference
            delta = ((1 / temperatures[j]) - (1 / temperatures[i])) * (
                lattices[i].calculate_energy() - lattices[j].calculate_energy()
            )

            if delta <= 0 or np.random.random() <= np.exp(-delta):
                lattices[i], lattices[j] = lattices[j], lattices[i]
            i += 2
        offset = 1 - offset
    # return the lattice with the lowest energy
    return min(lattices, key=lambda x: x.calculate_energy())
