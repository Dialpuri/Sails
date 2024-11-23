import argparse
from .__version__ import __version__
import importlib
from sails import morph, interface


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("-modelin", help="Path to input model", type=str, required=True)
    parser.add_argument(
        "-modelout",
        help="Path to output model",
        type=str,
        required=False,
        default="sails-model-out.cif",
    )
    parser.add_argument("-chain", help="Name of target chain", type=str, required=True)
    parser.add_argument(
        "-seqid", help="Name of target residue (protein)", type=int, required=True
    )
    parser.add_argument(
        "-wurcs", help="WURCS identifier for new glycan", type=str, required=True
    )

    return parser.parse_args()


def run():
    args = parse_args()

    sails_structure = interface.get_sails_structure(args.modelin)
    resource = importlib.resources.files("sails").joinpath("data")

    morphed_structure = morph(
        sails_structure, args.wurcs, args.chain, args.res, str(resource)
    )
    structure = interface.extract_sails_structure(morphed_structure)
    structure.make_mmcif_block().write_file(args.modelout)
