import argparse
from typing import Tuple
from sails import interface, n_glycosylate_from_objects, Dot, GlycoSite
import time
import gemmi
from pathlib import Path
import graphviz


def n_glycosylate(structure: gemmi.Structure | Path | str, mtz: gemmi.Mtz | Path | str, cycles: int, f: str, sigf: str) -> Tuple[
    gemmi.Structure, gemmi.Mtz]:
    sails_structure = get_sails_structure(structure)
    sails_mtz = get_sails_mtz(mtz, f, sigf)
    result = n_glycosylate_from_objects(sails_structure, sails_mtz, cycles)

    return interface.extract_sails_structure(result.structure), interface.extract_sails_mtz(result.mtz)


def get_sails_mtz(mtz: gemmi.Mtz | Path| str, f: str, sigf: str):
    if isinstance(mtz, gemmi.Mtz):
        sails_mtz = interface.extract_gemmi_mtz(mtz=mtz, column_names=[f, sigf])
    elif isinstance(mtz, Path) or isinstance(mtz, str):
        m = gemmi.read_mtz_file(str(mtz))
        sails_mtz = interface.extract_gemmi_mtz(mtz=m, column_names=[f, sigf])
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

    s, m = n_glycosylate(ip, im, 8)

    s.make_mmcif_block().write_file("sails-5fji.cif")
    # m.write_to_file("sails-5fji.mtz")

    # s = get_sails_structure(s)
    # d = Dot(s)
    # a = d.get_all_dotfiles()
    # for k, v in a.items():
    #     r = ip[k.model_idx][k.chain_idx][k.residue_idx]
    #     src = graphviz.Source(v)
    #     snfg = src.pipe(format='svg', encoding='utf-8')
    #     with open(f"testing/snfgs/{r.__str__()}.svg", "w") as f:
    #         f.write(snfg)

    t1 = time.time()
    print(f"Sails - Time Taken = {(t1 - t0)} seconds")


def run_cli():
    args = parse_args()
    t0 = time.time()

    f, sigf = args.colin.split(",")
    s, m = n_glycosylate(args.pdbin, args.mtzin, args.cycles, f, sigf)

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
    parser.add_argument("-colin", type=str, required=False, default="FP,SIGFP")
    parser.add_argument("-cycles", type=int, required=False, default=2)

    return parser.parse_args()
