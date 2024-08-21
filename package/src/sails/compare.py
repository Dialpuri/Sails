import argparse
import dataclasses
import importlib
from types import SimpleNamespace
import gemmi
import json
import functools


# Adapted from Paul Bond's Completeness Script in ModelCraft


def load_data_file(filename):
    with open(filename) as f:
        data = json.load(f)

    if not data:
        raise ValueError("Empty JSON file")
    residue_data = data["residues"]
    residue_data = {value["name"]: value for value in residue_data}
    return residue_data


def format_residue(chain: gemmi.Chain, residue: gemmi.Residue):
    return f"{chain.name}-{residue.name}-{residue.seqid}"


def match_residue(structure, search, radius, residue, atoms, same_type):
    return all(
        matching_atom(structure, search, radius, residue, name, same_type)
        for name in atoms
    )


def matching_atom(structure, search, radius, residue, name, same_type):
    for atom_alt in residue[name]:
        for mark in search.find_atoms(atom_alt.pos, "\0", radius=radius):
            cra = mark.to_cra(structure[0])
            if cra.atom.name == name:
                if not same_type or cra.residue.name == residue.name:
                    return True
    return False


@dataclasses.dataclass
class Statistics:
    sugars_total: int = 0
    sugars_built: int = 0


def main(args):
    resource = importlib.resources.files("sails").joinpath("data")
    data_file_path = resource / "data.json"
    data = load_data_file(data_file_path)
    structure = gemmi.read_structure(args.model)
    reference = gemmi.read_structure(args.reference)
    radius = args.radius
    search = gemmi.NeighborSearch(structure, max_radius=1).populate()

    matcher = functools.partial(match_residue, structure, search, radius)

    statistics = Statistics()

    for chain in reference[0]:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if residue.name not in data or info.is_amino_acid():
                continue

            ring_atoms = ["C1", "C2", "C3", "C4", "C5", "O5"]
            if any(name not in residue for name in ring_atoms):
                continue

            statistics.sugars_total += 1

            if matcher(residue, ring_atoms, True):
                statistics.sugars_built += 1

    return statistics


def compare(model_path: str, reference_path: str, radius: float):
    namespace = SimpleNamespace(
        model=model_path, reference=reference_path, radius=radius
    )
    stats = main(namespace)
    return stats


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("-model", type=str, required=True)
    parser.add_argument("-reference", type=str, required=True)
    parser.add_argument("-radius", type=float, required=False, default=1.5)

    args = parser.parse_args()

    stats = main(args)
    print(stats)


if __name__ == "__main__":
    run()
