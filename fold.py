import sys

from protein import Protein
from lattice import Lattice
from MCsearch import MCsearch
from REMCsearch import REMCsearch
from parser import parse_args


def main(args):
    args = parse_args(args)

    if args.protein:
        protein = Protein(args.protein)
    elif args.file:
        with open(args.file, "r") as f:
            protein = Protein(f.read())
    del args.protein, args.file

    lattice = Lattice(protein, args.initial_lattice)
    del args.initial_lattice

    # print(f"Initial lattice with energy of {lattice.calculate_energy()}")
    # lattice.draw_grid()

    sub_command = args.subparser_name
    del args.subparser_name

    print(args)
    if sub_command == "MC":
        final_lattice = MCsearch(**vars(args), lattice_input=lattice)
    elif sub_command == "REMC":
        final_lattice = REMCsearch(**vars(args), lattice_input=lattice)

    # print(f"Final lattice with energy of {lattice.calculate_energy()}")
    # final_lattice.draw_grid()


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
