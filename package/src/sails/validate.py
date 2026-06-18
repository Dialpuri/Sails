import argparse
import json
import time

from .__version__ import __version__
import importlib
from sails import validate, validate_site, interface

from .glycosylate import get_column_labels


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)

    parser.add_argument("--version", action="version", version=__version__)

    parent = argparse.ArgumentParser(add_help=False)
    group = parent.add_argument_group("Required arguments for all modes")
    group.add_argument("-v", action=argparse.BooleanOptionalAction, default=False)
    group.add_argument("--modelin", type=str, required=True)
    group.add_argument("--modelout", type=str, default="sails-validate.cif")
    group.add_argument("--logout", type=str, default="sails-validate.log")
    group.add_argument(
        "--threshold", type=float, default=0.8, help="RSCC Threshold to use for removal"
    )
    group.add_argument("--remove", action=argparse.BooleanOptionalAction, default=False)
    group.add_argument("--print", action=argparse.BooleanOptionalAction, default=False)

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
    xray_parser_group.add_argument("--chain", type=str, required=False)
    xray_parser_group.add_argument("--seqid", type=str, required=False)

    em_parser = subparsers.add_parser("em", parents=[parent], formatter_class=formatter)
    em_parser_group = em_parser.add_argument_group("Required arguments in EM mode")
    em_parser_group.add_argument("--mapin", type=str, required=True)
    em_parser_group.add_argument("--resolution", type=float, required=True)
    em_parser_group.add_argument(
        "--score", choices=["q", "rscc"], required=False, default="q"
    )

    return parser.parse_args()


def xray(args):
    sails_structure = interface.get_sails_structure(args.modelin)
    resource = importlib.resources.files("sails").joinpath("data")

    labels = get_column_labels(args.colin_fo, args.colin_fwt)
    sails_mtz = interface.get_sails_mtz(args.mtzin, *labels)

    if args.chain and args.seqid:
        result = validate_site(
            sails_structure,
            sails_mtz,
            args.chain,
            args.seqid,
            args.remove,
            args.threshold,
            str(resource),
        )
    else:
        result = validate(
            sails_structure, sails_mtz, args.remove, args.threshold, str(resource)
        )

    structure = interface.extract_sails_structure(result.structure)
    structure.make_mmcif_block().write_file(args.modelout)
    log = json.loads(result.log)

    if args.print:
        print(json.dumps(log, indent=4))

    with open(args.logout, "w") as f:
        json.dump(log, f, indent=4)


def em(args):
    sails_structure = interface.get_sails_structure(args.modelin)
    sails_grid = interface.get_sails_map(args.mapin)
    resource = importlib.resources.files("sails").joinpath("data")

    result = validate(
        sails_structure,
        sails_grid,
        args.resolution,
        args.remove,
        args.threshold,
        args.score == "q",
        str(resource),
    )

    structure = interface.extract_sails_structure(result.structure)
    structure.make_mmcif_block().write_file(args.modelout)
    log = json.loads(result.log)

    if args.print:
        print(json.dumps(log, indent=4))

    with open(args.logout, "w") as f:
        json.dump(log, f, indent=4)


def run():
    t0 = time.time()
    args = parse_args()

    if args.mode == "xray":
        xray(args)
    elif args.mode == "em":
        em(args)
    else:
        raise RuntimeError("Unknown mode")

    t1 = time.time()
    print(f"Sails Validate - Time Taken = {(t1 - t0)} seconds")
