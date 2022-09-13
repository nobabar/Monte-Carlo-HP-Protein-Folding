import argparse


class MyArgumentParser(argparse.ArgumentParser):
    def __innit__(self, *args, **kwargs):
        super().__innit__(self, *args, **kwargs)

    def print_help(self):
        # show main help
        help_text = self.format_help()

        # retrieve the subparsers
        subparsers_actions = [
            action
            for action in self._actions
            if isinstance(action, argparse._SubParsersAction)
        ]

        # show the complete help for the subparsers
        for subparsers_action in subparsers_actions:
            for choice, subparser in subparsers_action.choices.items():
                help_text += f"Subparser '{choice}'"
                help_text += subparser.format_help()

        print(help_text)


def parse_args(args):
    parser = MyArgumentParser()

    subparsers = parser.add_subparsers(help="sub-command help", dest="subparser_name")

    # parse the protein
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--protein", help="input protein sequence")
    group.add_argument(
        "-f", "--file", help="input file containing the protein sequence"
    )

    # create the parser for the Monte-Carlo command
    parser_MC = subparsers.add_parser("MC", help="Run the Monte Carlo algorithm")
    parser_MC.add_argument("-n", "--nsteps", type=int, default=1000)
    parser_MC.add_argument("-t", "--temperature", type=float, default=200.0)

    # create the parser for the Replica Exchange Monte-Carlo command
    parser_REMC = subparsers.add_parser(
        "REMC", help="Run the Replica Exchange Monte Carlo algorithm"
    )
    parser_REMC.add_argument("-n", "--nreplica", type=int, default=5)
    parser_REMC.add_argument("-e", "--energy-cutoff", type=int, default=-10)
    parser_REMC.add_argument("-m", "--max-steps", type=int, default=1000)
    parser_REMC.add_argument("-l", "--nlocal-steps", type=int, default=100)
    parser_REMC.add_argument("-tmin", "--temperature-min", type=float, default=100.0)
    parser_REMC.add_argument("-tmax", "--temperature-max", type=float, default=200.0)

    return parser.parse_args(args)


# args = parse_args("-p PHPPHPHPHPPHPPHPPHPHPPHPPHPHPPHP MC -n 1000 -t 200".split())
