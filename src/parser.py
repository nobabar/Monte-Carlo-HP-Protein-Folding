import argparse


class MyArgumentParser(argparse.ArgumentParser):
    """
    Parser class to override the help message
    """

    def __innit__(self, *args, **kwargs):
        super().__innit__(self, *args, **kwargs)

    def print_help(self, file=None):
        """
        Print main help message and subparsers help message together.
        Override the default help message.
        """
        # show main help
        help_text = self.format_help() + "\n"

        # retrieve the subparsers
        subparsers_actions = [
            action
            for action in self._actions
            if isinstance(action, argparse._SubParsersAction)
        ]

        # show the complete help for the subparsers
        for subparsers_action in subparsers_actions:
            for choice, subparser in subparsers_action.choices.items():
                help_text += f"Sub-command '{choice}'\n"
                help_text += subparser.format_help()

        print(help_text)


def parse_args(args):
    parser = MyArgumentParser()

    subparsers = parser.add_subparsers(
        help="sub-commands", dest="subparser_name")

    # parse the protein
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--protein", help="input protein sequence")
    group.add_argument(
        "-f", "--file", help="input file containing the protein sequence"
    )

    parser.add_argument("-i", "--initial-lattice", choices=["linear", "random"], default="linear",
                        help="initial lattice placement type, either in linear or using random walk")

    # create the parser for the Monte-Carlo command
    parser_MC = subparsers.add_parser(
        "MC", help="Run the Monte Carlo algorithm")
    parser_MC.add_argument("-n", "--n-steps", type=int, default=1000,
                           help="number of iterations in the search")
    parser_MC.add_argument("-t", "--temperature", type=float, default=200.0,
                           help="temperature of the search, higher temperatures will lead to more random movements")

    # create the parser for the Replica Exchange Monte-Carlo command
    parser_REMC = subparsers.add_parser(
        "REMC", help="Run the Replica Exchange Monte Carlo algorithm"
    )
    parser_REMC.add_argument("-n", "--n-replica", type=int,
                             default=5, help="number of replicas to use")
    parser_REMC.add_argument("-e", "--energy-cutoff", type=int, default=-10,
                             help="optimal energy to reach")
    parser_REMC.add_argument("-m", "--max-steps", type=int, default=1000,
                             help="maximum number of steps to perform if the energy cutoff is not reached")
    parser_REMC.add_argument("-l", "--local-steps", type=int, default=100,
                             help="number of steps to perform for each MC search")
    parser_REMC.add_argument("-tmin", "--temperature-min", type=float,
                             default=160.0, help="temperature of the first replica")
    parser_REMC.add_argument("-tmax", "--temperature-max", type=float,
                             default=220.0, help="temperature of the last replica")

    return parser.parse_args(args)
