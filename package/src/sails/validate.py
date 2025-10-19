import argparse
from .__version__ import __version__
import importlib
from sails import validate, interface

from .glycosylate import get_column_labels


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)

    parser.add_argument("--version", action="version", version=__version__)

    parent = argparse.ArgumentParser(add_help=False)
    group = parent.add_argument_group("Required arguments for all modes")
    group.add_argument("-v", action=argparse.BooleanOptionalAction, default=False)
    group.add_argument("--modelin", type=str, required=True)
    group.add_argument("--modelout", type=str, default="sails-validated.cif")
    group.add_argument("--remove", action=argparse.BooleanOptionalAction, default=False)

    formatter = argparse.ArgumentDefaultsHelpFormatter
    xray_parser = subparsers.add_parser(
        "xray", parents=[parent], formatter_class=formatter
    )
    xray_parser_group = xray_parser.add_argument_group(
        "Required arguments in X-ray mode"
    )
    xray_parser_group.add_argument("--mtzin", type=str, required=True)
    xray_parser_group.add_argument(
        "--colin-fo", type=str, required=False, default="FP,SIGFP"
    )
    xray_parser_group.add_argument(
        "--colin-fwt", type=str, required=False, default="FWT,PHWT"
    )

    em_parser = subparsers.add_parser("em", parents=[parent], formatter_class=formatter)
    em_parser_group = em_parser.add_argument_group("Required arguments in EM mode")
    em_parser_group.add_argument("--mapin", type=str, required=True)

    return parser.parse_args()


def run():
    args = parse_args()

    sails_structure = interface.get_sails_structure(args.modelin)
    resource = importlib.resources.files("sails").joinpath("data")

    labels = get_column_labels(args.colin_fo, args.colin_fwt)
    sails_mtz = interface.get_sails_mtz(args.mtzin, *labels)

    morphed_structure = validate(sails_structure, sails_mtz, args.remove, str(resource))
    structure = interface.extract_sails_structure(morphed_structure)
    structure.make_mmcif_block().write_file(args.modelout)
