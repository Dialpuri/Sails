import argparse

from .__version__ import __version__
import importlib
from sails import find_all_wurcs, find_wurcs, model_wurcs, interface
import json


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)

    parser.add_argument("--version", action="version", version=__version__)
    parent = argparse.ArgumentParser(add_help=False)
    group = parent.add_argument_group("Required arguments for all modes")
    # group.add_argument("-v", action=argparse.BooleanOptionalAction, default=False)
    group.add_argument("-modelin", type=str, required=True)

    formatter = argparse.ArgumentDefaultsHelpFormatter
    find_parser = subparsers.add_parser(
        "find", parents=[parent], formatter_class=formatter
    )
    find_parser.add_argument(
        "--all",
        help="Find WURCS identifiers for all glycans",
        action=argparse.BooleanOptionalAction,
        required=False,
    )
    find_parser.add_argument(
        "-chain", help="Name of target chain", type=str, required=False
    )
    find_parser.add_argument(
        "-seqid", help="Name of target residue (protein)", type=str, required=False
    )
    find_parser.add_argument(
        "-logout",
        help="Path to output json file",
        type=str,
        default="wurcs.json",
        required=False,
    )

    model_parser = subparsers.add_parser(
        "model", parents=[parent], formatter_class=formatter
    )
    model_parser.add_argument(
        "-modelout",
        help="Path to the output model",
        type=str,
        required=False,
        default="sails-model-out.cif",
    )
    model_parser.add_argument(
        "-wurcs",
        help="WURCS code to add to specified chain and residue",
        type=str,
        required=True,
    )
    model_parser.add_argument(
        "-chain", help="Name of target chain", type=str, required=True
    )
    model_parser.add_argument(
        "-seqid",
        help="Sequence ID of the root residue (protein)",
        type=int,
        required=True,
    )

    return parser.parse_args()


def run():
    args = parse_args()

    sails_structure = interface.get_sails_structure(args.modelin)
    resource = importlib.resources.files("sails").joinpath("data")

    if args.mode == "find":
        if not args.all and not args.chain:
            raise ValueError("Specify a chain and seqid, or the --all flag")
        if args.all:
            wurcs = find_all_wurcs(sails_structure, str(resource))
        else:
            wurcs = find_wurcs(sails_structure, args.chain, args.seqid, str(resource))
        with open(args.logout, "w", encoding="UTF-8") as json_file:
            json.dump(wurcs, json_file, indent=4)
    else:
        result = model_wurcs(
            sails_structure, args.wurcs, args.chain, args.seqid, str(resource)
        )
        structure = interface.extract_sails_structure(result)
        structure.make_mmcif_block().write_file(args.modelout)
