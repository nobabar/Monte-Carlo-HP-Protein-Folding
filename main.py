from protein import Protein
from lattice import Lattice
from MCsearch import MCsearch

protein = Protein("PHPPHPHPHPPHPPHPPHPHPPHPPHPHPPHP")
lattice = Lattice(protein, "linear")
print("Initial lattice:")
lattice.draw_grid()
lattice = MCsearch(1000, lattice, 200)
print("Final lattice:")
lattice.draw_grid()
