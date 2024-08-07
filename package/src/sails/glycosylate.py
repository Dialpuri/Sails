import argparse
import logging
from typing import Tuple, Callable, List
from sails import (interface, n_glycosylate_from_objects, c_glycosylate_from_objects, Dot, GlycoSite, __version__)
import time
import gemmi
from pathlib import Path
import json
import importlib.resources


def glycosylate(structure: gemmi.Structure | Path | str, mtz: gemmi.Mtz | Path | str, cycles: int, f: str, sigf: str,
                fwt: str, phwt: str, func: Callable = n_glycosylate_from_objects, verbose: bool = False) -> Tuple[
    gemmi.Structure, gemmi.Mtz, dict, dict]:
    """
    :param structure: The input structure file in gemmi.Structure, Path, or str format.
    :param mtz: The input MTZ file in gemmi.Mtz, Path, or str format.
    :param cycles: The number of cycles to perform glycosylation.
    :param f: The column label for the structure factor values.
    :param sigf: The column label for the structure factor uncertainties.
    :param fwt: The column label for the structure factor amplitude weights (potentially None).
    :param phwt: The column label for the structure factor phase weights (potentially None).
    :param func: The function to use for glycosylation. Default is n_glycosylate_from_objects.
    :param verbose: Flag specifying whether to print verbose output. Default is False.
    :return: A tuple containing the glycosylated structure in gemmi.Structure format,
             the glycosylated MTZ file in gemmi.Mtz format, the log as a string, and a dictionary of snfgs.
    """
    sails_structure = interface.get_sails_structure(structure)
    sails_mtz = interface.get_sails_mtz(mtz, f, sigf, fwt, phwt)
    resource = importlib.resources.files('sails').joinpath("data")
    result = func(sails_structure, sails_mtz, cycles, str(resource), verbose)

    return (interface.extract_sails_structure(result.structure), interface.extract_sails_mtz(result.mtz),
            json.loads(result.log), result.snfgs)


def get_column_labels(fo_columns: str, fwt_columns: str) -> List[str]:
    if ',' not in fo_columns:
        raise RuntimeError(f"Supplied colin-fo must be in form F,SIGF. Sails received {fo_columns}")

    fo_labels = [c for c in fo_columns.split(",") if c]
    fo_label_length = len(fo_labels)
    if fo_label_length != 2:
        raise RuntimeError(f"Too {'few' if fo_label_length < 2 else 'many'} column labels were provided")

    if fwt_columns == "":
        return fo_labels + [None, None]

    if ',' not in fwt_columns:
        raise RuntimeError(f"Supplied colin-fwt must be in form FWT,PHWT. Sails received {fo_columns}")

    fwt_labels = [c for c in fwt_columns.split(",") if c]
    fwt_label_length = len(fwt_labels)
    if fwt_label_length != 2:
        raise RuntimeError(f"Too {'few' if fwt_label_length < 2 else 'many'} column labels were provided")

    return fo_labels + fwt_labels


def save_log(log: dict, args: argparse.Namespace) -> None:
    arguments = vars(args)
    log['arguments'] = arguments
    with open(args.logout, "w") as f:
        json.dump(log, f, indent=4)


def save_snfgs(snfgs: dict, snfg_path: Path):
    output_path = snfg_path
    if snfg_path.suffixes:
        logging.warning("Supplied SNFG path is a file, not a directory. A new directory will be made.")
        output_path = snfg_path.parent / "sails-snfgs"

    output_path.mkdir(exist_ok=True)

    for cycle_no, snfg_strings in snfgs.items():
        cycle = output_path / f"cycle-{cycle_no}"
        cycle.mkdir(exist_ok=True)
        for site, snfg_string in snfg_strings.items():
            path = cycle / f"{site}.svg"
            with open(path, "w") as f:
                f.write(snfg_string)


def run_cli():
    args = parse_args()
    t0 = time.time()

    labels = get_column_labels(args.colin_fo, args.colin_fwt)

    func = c_glycosylate_from_objects if args.cglycan else n_glycosylate_from_objects
    cycles = 1 if args.cglycan else args.cycles
    structure, mtz, log, snfgs = glycosylate(args.pdbin, args.mtzin, cycles, *labels, func, args.v)

    if args.snfgout:
        save_snfgs(snfgs, Path(args.snfgout))

    structure.make_mmcif_block().write_file(args.pdbout)
    mtz.write_to_file(args.mtzout)

    save_log(log, args)

    t1 = time.time()
    print(f"Sails - Time Taken = {(t1 - t0)} seconds")


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-pdbin", type=str, required=True)
    parser.add_argument("-mtzin", type=str, required=True)
    parser.add_argument("-pdbout", type=str, required=False, default="sails-model-out.cif")
    parser.add_argument("-mtzout", type=str, required=False, default="sails-refln-out.mtz")
    parser.add_argument("-logout", type=str, default="sails-log.json")
    parser.add_argument("-snfgout", type=str)
    parser.add_argument("-colin-fo", type=str, required=False, default="FP,SIGFP")
    parser.add_argument("-colin-fwt", type=str, required=False, default="")
    parser.add_argument("-cycles", type=int, required=False, default=2)
    parser.add_argument("-nglycan", action=argparse.BooleanOptionalAction)
    parser.add_argument("-cglycan", action=argparse.BooleanOptionalAction)
    parser.add_argument("-v", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--version", action="version", version=__version__)

    return parser.parse_args()
