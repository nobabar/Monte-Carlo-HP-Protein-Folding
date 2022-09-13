#!/usr/bin/env python

import sys

from protein import Protein
from lattice import Lattice
from MCsearch import MCsearch
from REMCsearch import REMCsearch
from parser import parse_args


def main(args):
    args = parse_args(args)

    print(args)

    # protein = Protein("PHPPHPHPHPPHPPHPPHPHPPHPPHPHPPHP")
    # lattice = Lattice(protein, "random")

    # print(f"Initial lattice with energy of {lattice.calculate_energy()}")
    # lattice.draw_grid()

    # # lattice = MCsearch(5000, lattice, 500)
    # lattice = REMCsearch(5, -5, 500, 100, 160, 220, lattice)

    # print(f"Final lattice with energy of {lattice.calculate_energy()}")
    # lattice.draw_grid()


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
