from collections import defaultdict
from pathlib import Path
from typing import List, Tuple
import gemmi
import argparse
import json


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


def run():
    """
    Parse command-line arguments, read PDB model, find glycosylation sites,
    and write the results to an output file in JSON format.

    :return: None
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-modelin", required=True, type=str, help="Path to a model in PDB or CIF format"
    )
    parser.add_argument(
        "-logout",
        required=False,
        default="sites.json",
        type=str,
        help="Path to output file",
    )

    args = parser.parse_args()

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
