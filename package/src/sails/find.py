import importlib
import time
from argparse import ArgumentError
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple
import gemmi
import argparse
import json
from sails import identify_predicted_sites, GlycoSite
from .interface import get_sails_structure, get_sails_map
from .glycosylate import read_prediction_dir, save_log
from .prediction.model import ModelType
from .prediction.predict import predict_map


def find_n_glycosylation_sites(structure: gemmi.Structure):
    """
    :param structure: A gemmi.Structure object representing the structure to search for N-glycosylation sites.
    :return: A list of tuples, where each tuple represents the position of an N-glycosylation site in the structure.
     Each tuple contains the following elements:
        - mi: The model index of the residue containing the N-glycosylation site.
        - ci: The chain index of the residue containing the N-glycosylation site.
        - ri: The residue index of the N-glycosylation site within the chain.
    """
    sites = []
    for mi, m in enumerate(structure):
        for ci, c in enumerate(m):
            if len(c) < 3:
                continue
            for ri in range(len(c) - 2):
                first = gemmi.find_tabulated_residue(c[ri].name).one_letter_code
                if first != "N":
                    continue

                third = gemmi.find_tabulated_residue(c[ri + 2].name).one_letter_code
                if third != "S" and third != "T":
                    continue

                second = gemmi.find_tabulated_residue(c[ri + 1].name).one_letter_code
                if second == "P":
                    continue

                sites.append((mi, ci, ri))
    return sites


def find_c_glycosylation_sites(structure: gemmi.Structure):
    """
    Finds C-glycosylation sites in a given structure.

    :param structure: a gemmi.Structure object representing the structure to search in
    :return: a list of tuples representing the indices of the C-glycosylation sites
    Each tuple contains the following elements:
        - mi: The model index of the residue containing the N-glycosylation site.
        - ci: The chain index of the residue containing the N-glycosylation site.
        - ri: The residue index of the N-glycosylation site within the chain.
    """
    sites = []
    for mi, m in enumerate(structure):
        for ci, c in enumerate(m):
            if len(c) < 4:
                continue
            for ri in range(len(c) - 3):
                first = gemmi.find_tabulated_residue(c[ri].name).one_letter_code
                if first != "W":
                    continue

                fourth = gemmi.find_tabulated_residue(c[ri + 3].name).one_letter_code
                if fourth != "W":
                    continue

                sites.append((mi, ci, ri))
                sites.append((mi, ci, ri + 3))

    return sites


def format_sites(
    sites: List[Tuple[int, int, int]], structure: gemmi.Structure
) -> List[dict]:
    """Format the given sites data.

    :param sites: A list of tuples representing the site data.
    :param structure: The gemmi.Structure object containing the site data.
    :return: A list of dictionaries containing the formatted site data.
    """
    d = []
    for site in sites:
        mi, ci, ri = site
        c = structure[mi][ci]
        r = c[ri]
        entry = {
            "chainId": c.name,
            "residueSeqId": r.seqid.__str__(),
            "residueName": r.name,
            "modelIndex": mi,
            "chainIndex": ci,
            "residueIndex": ri,
        }
        d.append(entry)
    return d


def sequence_find(args: argparse.Namespace):
    pdb_path = Path(args.modelin)
    if not pdb_path.exists():
        raise FileNotFoundError("Could not find specified file")

    structure = gemmi.read_structure(args.modelin)
    data = defaultdict(list)

    n_glycosylation_sites = find_n_glycosylation_sites(structure)
    if n_glycosylation_sites:
        data["n-glycosylation"] = format_sites(n_glycosylation_sites, structure)

    c_glycosylation_sites = find_c_glycosylation_sites(structure)
    if c_glycosylation_sites:
        data["c-glycosylation"] = format_sites(c_glycosylation_sites, structure)

    with open(args.logout, "w") as f:
        json.dump(data, f, indent=4)


def convert_residue_name_to_type(residue_name: str) -> str:
    n_glycans = ["ASN"]
    o_glycans = ["SER", "THR"]
    c_glycans = ["TRP"]

    if residue_name in n_glycans:
        return "n-glycan"
    elif residue_name in o_glycans:
        return "o-glycan"
    elif residue_name in c_glycans:
        return "c-glycan"
    return "x-glycan"


def convert_glycosites_to_log(
    glycosites: List[GlycoSite], structure: gemmi.Structure | Path | str
):
    if isinstance(structure, str) or isinstance(structure, Path):
        structure = gemmi.read_structure(str(structure))

    keys = defaultdict(list)
    for glycosite in glycosites:
        model = structure[glycosite.model_idx]
        chain = model[glycosite.chain_idx]
        residue = chain[glycosite.residue_idx]
        key = f"{chain.name}-{residue.name}-{residue.seqid.num}"
        keys[convert_residue_name_to_type(residue.name)].append(key)

    return keys


def get_amplitude_phase(args):
    if "," not in args.colin_fwt:
        raise ArgumentError("FWT column should be comma separated")
    return args.colin_fwt.split(",")


def xray(args):
    sails_structure = get_sails_structure(args.modelin)
    resource = importlib.resources.files("sails").joinpath("data")
    model = ModelType[args.modeltype]
    if args.preddirin:
        predictions = read_prediction_dir(args.preddirin, model)
    else:
        amplitude, phase = get_amplitude_phase(args)
        predictions = predict_map(
            model.name,
            args.mtzin,
            "output",
            nthreads=8,
            amplitude=amplitude,
            phase=phase,
            save_map=True,
        )

    if model == ModelType.binary:
        glycan_predicted_map = predictions
        sails_grid = get_sails_map(glycan_predicted_map)
        result = identify_predicted_sites(sails_structure, sails_grid, str(resource))
    else:
        glycan_predicted_map, protein_predicted_map = predictions
        sails_glycan_grid = get_sails_map(glycan_predicted_map)
        sails_protein_grid = get_sails_map(protein_predicted_map)
        searchtype = args.searchtype
        result = identify_predicted_sites(
            sails_structure,
            sails_glycan_grid,
            sails_protein_grid,
            searchtype == "glycan",
            str(resource),
        )

    log = convert_glycosites_to_log(result, args.modelin)
    save_log(log, args)


def em(args):
    sails_structure = get_sails_structure(args.modelin)
    resource = importlib.resources.files("sails").joinpath("data")
    model = ModelType[args.modeltype]

    if args.preddirin:
        predictions = read_prediction_dir(args.preddirin, model)
    else:
        predictions = predict_map(
            model.name,
            args.mapin,
            "output",
            nthreads=8,
            save_map=True,
        )

    if model == ModelType.binary:
        glycan_predicted_map = predictions
        sails_grid = get_sails_map(glycan_predicted_map)
        result = identify_predicted_sites(sails_structure, sails_grid, str(resource))
    else:
        glycan_predicted_map, protein_predicted_map = predictions
        sails_glycan_grid = get_sails_map(glycan_predicted_map)
        sails_protein_grid = get_sails_map(protein_predicted_map)
        searchtype = args.searchtype
        result = identify_predicted_sites(
            sails_structure,
            sails_glycan_grid,
            sails_protein_grid,
            searchtype == "glycan",
            str(resource),
        )

    log = convert_glycosites_to_log(result, args.modelin)
    save_log(log, args)


def density_find(args: argparse.Namespace):
    t0 = time.time()

    if args.source == "xray":
        xray(args)
    elif args.source == "em":
        em(args)
    else:
        raise RuntimeError("Unknown mode")

    t1 = time.time()
    print(f"Sails Density Identification - Time Taken = {(t1 - t0)} seconds")


def run():
    """
    Parse command-line arguments, read PDB model, find glycosylation sites,
    and write the results to an output file in JSON format.

    :return: None
    """

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest="mode", required=True)

    seq_parser = subparsers.add_parser("seq")
    seq_parser.add_argument(
        "--modelin",
        required=True,
        type=str,
        help="Path to a model in PDB or CIF format",
    )
    seq_parser.add_argument(
        "--logout",
        required=False,
        default="sites.json",
        type=str,
        help="Path to output file",
    )

    density_parser = subparsers.add_parser("density")
    density_subparser = density_parser.add_subparsers(dest="source", required=True)
    xray_parser = density_subparser.add_parser("xray")
    xray_parser.add_argument(
        "--mtzin", required=True, type=str, help="Path to mtz file"
    )
    xray_parser.add_argument(
        "--modelin",
        required=True,
        type=str,
        help="Path to a model in PDB or CIF format",
    )
    xray_parser.add_argument(
        "--preddirin",
        required=False,
        type=str,
        help="Path to a model in PDB or CIF format",
    )
    xray_parser.add_argument(
        "--logout",
        required=False,
        default="sites.json",
        type=str,
        help="Path to output file",
    )
    xray_parser.add_argument(
        "--modeltype",
        required=True,
        choices=[type.name for type in ModelType],
        help="Binary or Multiclass model",
    )
    xray_parser.add_argument(
        "--searchtype",
        required=True,
        choices=["protein", "glycan"],
        help="Search for protein or glycan, only used if modeltype is multiclass",
    )
    xray_parser.add_argument("--colin-fo", type=str, required=False, default="FP,SIGFP")
    xray_parser.add_argument(
        "--colin-fwt", type=str, required=False, default="FWT,PHWT"
    )

    em_parser = density_subparser.add_parser("em")
    em_parser.add_argument("--mapin", type=str, required=True)
    em_parser.add_argument(
        "--modelin",
        required=True,
        type=str,
        help="Path to a model in PDB or CIF format",
    )
    em_parser.add_argument(
        "--logout",
        required=False,
        default="sites.json",
        type=str,
        help="Path to output file",
    )
    em_parser.add_argument(
        "--preddirin",
        required=False,
        type=str,
        help="Path to a model in PDB or CIF format",
    )
    em_parser.add_argument(
        "--modeltype",
        required=True,
        choices=[type.name for type in ModelType],
        help="Binary or Multiclass model",
    )
    em_parser.add_argument(
        "--searchtype",
        required=True,
        choices=["protein", "glycan"],
        help="Search for protein or glycan, only used if modeltype is multiclass",
    )

    args = parser.parse_args()
    if args.mode == "seq":
        sequence_find(args)
    elif args.mode == "density":
        density_find(args)
