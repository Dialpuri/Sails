import argparse
from typing import Tuple, Callable
from sails import interface, n_glycosylate_from_objects, c_glycosylate_from_objects, Dot, GlycoSite
import time
import gemmi
from pathlib import Path
import json
import graphviz
import importlib.resources


def glycosylate(structure: gemmi.Structure | Path | str, mtz: gemmi.Mtz | Path | str, cycles: int, f: str, sigf: str,
                func: Callable = n_glycosylate_from_objects, verbose: bool = False) -> Tuple[
    gemmi.Structure, gemmi.Mtz, str]:
    """
    :param structure: The input structure file in gemmi.Structure, Path, or str format.
    :param mtz: The input MTZ file in gemmi.Mtz, Path, or str format.
    :param cycles: The number of cycles to perform glycosylation.
    :param f: The column label for the structure factor values.
    :param sigf: The column label for the structure factor uncertainties.
    :param func: The function to use for glycosylation. Default is n_glycosylate_from_objects.
    :param verbose: Flag specifying whether to print verbose output. Default is False.
    :return: A tuple containing the glycosylated structure in gemmi.Structure format,
             the glycosylated MTZ file in gemmi.Mtz format, and the log as a string.
    """
    sails_structure = get_sails_structure(structure)
    sails_mtz = get_sails_mtz(mtz, f, sigf)
    resource = importlib.resources.files('sails').joinpath("data")
    result = func(sails_structure, sails_mtz, cycles, str(resource), verbose)

    return (interface.extract_sails_structure(result.structure), interface.extract_sails_mtz(result.mtz),
            json.loads(result.log))


def read_sf_cif(mtz: Path):
    """
    :param mtz: Path to the CIF file containing structure factor data.
    :return: MTZ file containing the converted structure factor data.

    """
    doc = gemmi.cif.read(str(mtz))
    rblocks = gemmi.as_refln_blocks(doc)
    if not rblocks:
        raise RuntimeError("SF CIF supplied has no reflections")
    cif2mtz = gemmi.CifToMtz()
    return cif2mtz.convert_block_to_mtz(rblocks[0])


def get_sails_mtz(mtz: gemmi.Mtz | Path | str, f: str, sigf: str):
    """
    :param mtz: Path to an MTZ file, an instance of gemmi.Mtz, or a string representing the path to an MTZ file.
    :param f: Column name of the F values in the MTZ file.
    :param sigf: Column name of the SIGF values in the MTZ file.
    :return: A gemmi.Mtz object containing only the specified F and SIGF columns.

    This method extracts the specified F and SIGF columns from an MTZ file and returns a new MTZ object
    containing only those columns. It supports three types of input for the 'mtz' parameter:
    1) A string or Path object representing the path to the MTZ file.
    2) An instance of the gemmi.Mtz class.
    3) A gemmi.Mtz object.

    If 'mtz' is an instance of the gemmi.Mtz class, the method directly extracts the specified columns
    using the provided 'f' and 'sigf' column names.

    If 'mtz' is a string or Path object representing the path to an MTZ file, the method first checks
    the file extension to determine the format of the file. If the file has a CIF or ENT extension,
    it is assumed to be a CIF file and the 'read_sf_cif' function is used to read the data. Otherwise,
    the 'gemmi.read_mtz_file' function is used to read the MTZ file. The extracted columns are then
    obtained from the resulting MTZ object using the specified 'f' and 'sigf' column names.

    If 'mtz' is neither an instance of gemmi.Mtz, nor a string/Path object, a RuntimeError is raised.

    The method returns a new gemmi.Mtz object containing only the specified F and SIGF columns.
    """
    if isinstance(mtz, gemmi.Mtz):
        sails_mtz = interface.extract_gemmi_mtz(mtz=mtz, column_names=[f, sigf])
    elif isinstance(mtz, Path) or isinstance(mtz, str):
        mtz = Path(mtz)
        if ".cif" in mtz.suffixes or ".ent" in mtz.suffixes:
            m = read_sf_cif(mtz)
        else:
            m = gemmi.read_mtz_file(str(mtz))
        sails_mtz = interface.extract_gemmi_mtz(mtz=m, column_names=[f, sigf])
    else:
        raise RuntimeError("Unknown object passed to second argument of n_glycosylate function, allowed types are"
                           "gemmi.Mtz, Path, and str")
    return sails_mtz


def get_sails_structure(structure):
    """
    Retrieves the Sails structure from the provided input.

    :param structure: The input structure. It can be either of type gemmi.Structure, or a file path as type Path
                      or str.
    :return: The Sails structure object.
    :raises RuntimeError: If an unknown object is passed as the first argument. Only gemmi.Structure, Path, and str
                          are allowed types.
    """
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

    s, m, l = glycosylate(ip, im, 1, "FP", "SIGFP")

    with open("sails-5fji-log.json", "w") as f:
        json.dump(l, f, indent=4)

    s.make_mmcif_block().write_file("sails-5fji-test.cif")
    m.write_to_file("sails-5fji-test.mtz")

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
    if args.cglycan:
        s, m = glycosylate(args.pdbin, args.mtzin, args.cycles, f, sigf, c_glycosylate_from_objects, args.v)
    else:
        s, m = glycosylate(args.pdbin, args.mtzin, args.cycles, f, sigf, n_glycosylate_from_objects, args.v)

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
    parser.add_argument("-nglycan", action=argparse.BooleanOptionalAction)
    parser.add_argument("-cglycan", action=argparse.BooleanOptionalAction)
    parser.add_argument("-v", action=argparse.BooleanOptionalAction, default=False)

    return parser.parse_args()
