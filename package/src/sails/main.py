import argparse
from typing import Tuple
from sails import interface, n_glycosylate_from_objects, Dot, GlycoSite
import time
import gemmi
from pathlib import Path
import graphviz


def n_glycosylate(structure: gemmi.Structure | Path | str, mtz: gemmi.Mtz | Path | str, cycles: int) -> Tuple[
    gemmi.Structure, gemmi.Mtz]:
    sails_structure = get_sails_structure(structure)
    sails_mtz = get_sails_mtz(mtz)
    result = n_glycosylate_from_objects(sails_structure, sails_mtz, cycles)

    return interface.extract_sails_structure(result.structure), interface.extract_sails_mtz(result.mtz)


def get_sails_mtz(mtz):
    if isinstance(mtz, gemmi.Mtz):
        sails_mtz = interface.extract_gemmi_mtz(mtz=mtz)
    elif isinstance(mtz, Path) or isinstance(mtz, str):
        m = gemmi.read_mtz_file(str(mtz))
        sails_mtz = interface.extract_gemmi_mtz(mtz=m)
    else:
        raise RuntimeError("Unknown object passed to second argument of n_glycosylate function, allowed types are"
                           "gemmi.Mtz, Path, and str")
    return sails_mtz


def get_sails_structure(structure):
    if isinstance(structure, gemmi.Structure):
        sails_structure = interface.extract_gemmi_structure(structure=structure)
    elif isinstance(structure, Path) or isinstance(structure, str):
        s = gemmi.read_structure(str(structure))
        sails_structure = interface.extract_gemmi_structure(structure=s)
    else:
        raise RuntimeError("Unknown object passed to first argument of n_glycosylate function, allowed types are"
                           "gemmi.Structure, Path, and str")
    return sails_structure


def run_python():
    print("Running Sails")
    t0 = time.time()
    ip = gemmi.read_structure("package/models/5fji/5fji_deglycosylated.pdb")
    im = gemmi.read_mtz_file("package/models/5fji/5fji.mtz")
    s, m = n_glycosylate(ip, im, 3)

    s.make_mmcif_block().write_file("sails-5fji.cif")
    m.write_to_file("sails-5fji.mtz")

    # d = Dot(get_sails_structure(ip))
    # g = GlycoSite(0, 0, 40)
    # string = d.get_dotfile(g)
    # src = graphviz.Source(string)
    # snfg = src.pipe(format='svg', encoding='utf-8')
    # with open("output.svg", "w") as f:
    #     f.write(snfg)

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

    parser.add_argument("-pdbin", type=str, required=True)
    parser.add_argument("-mtzin", type=str, required=True)
    parser.add_argument("-pdbout", type=str, required=False)
    parser.add_argument("-mtzout", type=str, required=False)
    parser.add_argument("-cycles", type=int, required=False, default=2)

    return parser.parse_args()
