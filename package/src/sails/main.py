import argparse
from typing import Tuple
from sails import interface, n_glycosylate_from_objects
import time
import gemmi
from pathlib import Path


def n_glycosylate(structure: gemmi.Structure | Path | str, mtz: gemmi.Mtz | Path | str, cycles: int) -> Tuple[
    gemmi.Structure, gemmi.Mtz]:
    if isinstance(structure, gemmi.Structure):
        sails_structure = interface.extract_gemmi_structure(structure=structure)
    elif isinstance(structure, Path) or isinstance(structure, str):
        s = gemmi.read_structure(str(structure))
        sails_structure = interface.extract_gemmi_structure(structure=s)
    else:
        raise RuntimeError("Unknown object passed to first argument of n_glycosylate function, allowed types are"
                           "gemmi.Structure, Path, and str")

    if isinstance(mtz, gemmi.Mtz):
        sails_mtz = interface.extract_gemmi_mtz(mtz=mtz)
    elif isinstance(mtz, Path) or isinstance(mtz, str):
        m = gemmi.read_mtz_file(str(mtz))
        sails_mtz = interface.extract_gemmi_mtz(mtz=m)
    else:
        raise RuntimeError("Unknown object passed to second argument of n_glycosylate function, allowed types are"
                           "gemmi.Mtz, Path, and str")

    result = n_glycosylate_from_objects(sails_structure, sails_mtz, cycles)

    return interface.extract_sails_structure(result.structure), interface.extract_sails_mtz(result.mtz)


def run_python():
    print("Running Sails")
    t0 = time.time()
    ip = gemmi.read_structure("package/models/5fji/5fji_deglycosylated.pdb")
    im = gemmi.read_mtz_file("package/models/5fji/5fji.mtz")
    s, m = n_glycosylate(ip, im, 7)
    #
    s.make_mmcif_block().write_file("sails-5fji.cif")
    m.write_to_file("sails-5fji.mtz")

    t1 = time.time()
    print(f"Sails - Time Taken = {(t1 - t0)} seconds")


def run_cli():
    args = parse_args()
    t0 = time.time()

    s, m = n_glycosylate(args.pdbin, args.mtzin, args.cycles)

    s.make_mmcif_block().write_file(args.pdbout)
    m.write_to_file(args.mtzout)

    t1 = time.time()
    print(f"Sails - Time Taken = {(t1 - t0)} seconds")


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--pdbin", type=str, required=True)
    parser.add_argument("-m", "--mtzin", type=str, required=True)
    parser.add_argument("--pdbout", type=str, required=False)
    parser.add_argument("--mtzout", type=str, required=False)
    parser.add_argument("--cycles", type=int, required=False, default=2)

    return parser.parse_args()
