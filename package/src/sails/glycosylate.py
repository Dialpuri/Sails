import argparse
import enum
import importlib.resources
import json
import logging
import time
from pathlib import Path
from typing import Tuple, List

import gemmi
from sails import interface, n_glycosylate, c_glycosylate, o_mannosylate, __version__


class Type(enum.IntEnum):
    n_glycosylate = 1
    c_glycosylate = 2
    o_mannosylate = 3

    def __str__(self):
        return self.name

    @staticmethod
    def from_string(s: str):
        try:
            return Type[s]
        except KeyError:
            raise ValueError("Invalid Type")


def map_type_to_function(type: Type):
    if type == Type.n_glycosylate:
        return n_glycosylate

    if type == Type.c_glycosylate:
        return c_glycosylate

    if type == Type.o_mannosylate:
        return o_mannosylate

    raise TypeError("Type not found")


def glycosylate_xtal(
    structure: gemmi.Structure | Path | str,
    mtz: gemmi.Mtz | Path | str,
    cycles: int,
    f: str,
    sigf: str,
    fwt: str,
    phwt: str,
    type: Type = Type.n_glycosylate,
    verbose: bool = False,
) -> Tuple[gemmi.Structure, gemmi.Mtz, dict, dict]:
    """
    :param structure: The input structure file in gemmi.Structure, Path, or str format.
    :param mtz: The input MTZ file in gemmi.Mtz, Path, or str format.
    :param cycles: The number of cycles to perform glycosylation.
    :param f: The column label for the structure factor values.
    :param sigf: The column label for the structure factor uncertainties.
    :param fwt: The column label for the structure factor amplitude weights (potentially None).
    :param phwt: The column label for the structure factor phase weights (potentially None).
    :param type: The type of glycosylation to perform. Default is n_glycosylate.
    :param verbose: Flag specifying whether to print verbose output. Default is False.
    :return: A tuple containing the glycosylated structure in gemmi.Structure format,
             the glycosylated MTZ file in gemmi.Mtz format, the log as a string, and a dictionary of snfgs.
    """
    sails_structure = interface.get_sails_structure(structure)
    sails_mtz = interface.get_sails_mtz(mtz, f, sigf, fwt, phwt)
    resource = importlib.resources.files("sails").joinpath("data")

    func = map_type_to_function(type)
    result = func(sails_structure, sails_mtz, cycles, str(resource), verbose)

    return (
        interface.extract_sails_structure(result.structure),
        interface.extract_sails_mtz(result.mtz),
        json.loads(result.log),
        result.snfgs,
    )


def glycosylate_em(
    structure: gemmi.Structure | Path | str,
    map: gemmi.Ccp4Map | gemmi.FloatGrid | Path | str,
    cycles: int,
    type: Type = Type.n_glycosylate,
    verbose: bool = False,
) -> Tuple[gemmi.Structure, dict, dict]:
    sails_structure = interface.get_sails_structure(structure)
    sails_grid = interface.get_sails_map(map)
    resource = importlib.resources.files("sails").joinpath("data")

    func = map_type_to_function(type)
    result = func(sails_structure, sails_grid, cycles, str(resource), verbose)

    return (
        interface.extract_sails_structure(result.structure),
        json.loads(result.log),
        result.snfgs,
    )


def get_column_labels(fo_columns: str, fwt_columns: str) -> List[str]:
    if "," not in fo_columns:
        raise RuntimeError(
            f"Supplied colin-fo must be in form F,SIGF. Sails received {fo_columns}"
        )

    fo_labels = [c for c in fo_columns.split(",") if c]
    fo_label_length = len(fo_labels)
    if fo_label_length != 2:
        raise RuntimeError(
            f"Too {'few' if fo_label_length < 2 else 'many'} column labels were provided"
        )

    if fwt_columns == "":
        return fo_labels + [None, None]

    if "," not in fwt_columns:
        raise RuntimeError(
            f"Supplied colin-fwt must be in form FWT,PHWT. Sails received {fo_columns}"
        )

    fwt_labels = [c for c in fwt_columns.split(",") if c]
    fwt_label_length = len(fwt_labels)
    if fwt_label_length != 2:
        raise RuntimeError(
            f"Too {'few' if fwt_label_length < 2 else 'many'} column labels were provided"
        )

    return fo_labels + fwt_labels


def save_log(log: dict, args: argparse.Namespace) -> None:
    arguments = vars(args)
    log["arguments"] = arguments
    with open(args.logout, "w") as f:
        json.dump(log, f, indent=4)


def save_snfgs(snfgs: dict, snfg_path: Path):
    output_path = snfg_path
    if snfg_path.suffixes:
        logging.warning(
            "Supplied SNFG path is a file, not a directory. A new directory will be made."
        )
        output_path = snfg_path.parent / "sails-snfgs"

    output_path.mkdir(exist_ok=True)

    for cycle_no, snfg_strings in snfgs.items():
        cycle = output_path / f"cycle-{cycle_no}"
        cycle.mkdir(exist_ok=True)
        for site, snfg_string in snfg_strings.items():
            path = cycle / f"{site}.svg"
            with open(path, "w") as f:
                f.write(snfg_string)


def xray(args):
    labels = get_column_labels(args.colin_fo, args.colin_fwt)

    cycles = args.cycles if args.type == Type.n_glycosylate else 1
    structure, mtz, log, snfgs = glycosylate_xtal(
        args.modelin, args.mtzin, cycles, *labels, args.type, args.v
    )

    if args.snfgout:
        save_snfgs(snfgs, Path(args.snfgout))

    structure.make_mmcif_block().write_file(args.modelout)
    mtz.write_to_file(args.mtzout)

    save_log(log, args)


def em(args):
    cycles = args.cycles if args.type == Type.n_glycosylate else 1
    structure, log, snfgs = glycosylate_em(
        args.modelin, args.mapin, cycles, args.type, args.v
    )
    structure.make_mmcif_block().write_file(args.modelout)
    save_log(log, args)


def run_cli():
    t0 = time.time()
    args = parse_args()

    if args.mode == "xray":
        xray(args)
    elif args.mode == "em":
        em(args)
    else:
        raise RuntimeError("Unknown mode")

    t1 = time.time()
    print(f"Sails - Time Taken = {(t1 - t0)} seconds")


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)

    parser.add_argument("--version", action="version", version=__version__)

    parent = argparse.ArgumentParser(add_help=False)
    group = parent.add_argument_group("Required arguments for all modes")
    group.add_argument("-v", action=argparse.BooleanOptionalAction, default=False)
    group.add_argument("-modelin", type=str, required=True)
    group.add_argument(
        "-modelout", type=str, required=False, default="sails-model-out.cif"
    )
    group.add_argument("-logout", type=str, default="sails-log.json")
    group.add_argument("-snfgout", type=str)
    group.add_argument("-cycles", type=int, required=False, default=2)
    group.add_argument(
        "-type", type=Type.from_string, choices=list(Type), default=Type.n_glycosylate
    )

    formatter = argparse.ArgumentDefaultsHelpFormatter
    xray_parser = subparsers.add_parser(
        "xray", parents=[parent], formatter_class=formatter
    )
    xray_parser_group = xray_parser.add_argument_group(
        "Required arguments in X-ray mode"
    )
    xray_parser_group.add_argument("-mtzin", type=str, required=True)
    xray_parser_group.add_argument(
        "-mtzout", type=str, required=False, default="sails-refln-out.mtz"
    )
    xray_parser_group.add_argument(
        "-colin-fo", type=str, required=False, default="FP,SIGFP"
    )
    xray_parser_group.add_argument("-colin-fwt", type=str, required=False, default="")

    em_parser = subparsers.add_parser("em", parents=[parent], formatter_class=formatter)
    em_parser_group = em_parser.add_argument_group("Required arguments in EM mode")
    em_parser_group.add_argument("-mapin", type=str, required=True)

    return parser.parse_args()
