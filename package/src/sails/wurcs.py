import argparse
from .__version__ import __version__
import importlib
from sails import wurcs, interface


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("-modelin", help="Path to input model", type=str, required=True)
    parser.add_argument("-chain", help="Name of target chain", type=str, required=True)
    parser.add_argument(
        "-res", help="Name of target residue (protein)", type=str, required=True
    )

    return parser.parse_args()


def run():
    args = parse_args()

    sails_structure = interface.get_sails_structure(args.modelin)
    resource = importlib.resources.files("sails").joinpath("data")

    wurcs(sails_structure, args.chain, args.res, str(resource))
