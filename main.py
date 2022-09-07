from protein import Protein
from lattice import Lattice
from MCsearch import MCsearch

protein = Protein("PHPPHPHPHPPHPPHPPHPHPPHPPHPHPPHP")
lattice = Lattice(protein, "linear")

print(f"Initial lattice with energy of {lattice.calculate_energy()}")
lattice.draw_grid()

lattice = MCsearch(1000, lattice, 300)

print(f"Final lattice with energy of {lattice.calculate_energy()}")
lattice.draw_grid()
