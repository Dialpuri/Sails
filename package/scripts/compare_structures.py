import argparse
import gemmi
import json
import numpy as np


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


def main(args):
    data = load_data_file("package/data/data.json")
    structure = gemmi.read_structure(args.model)
    reference = gemmi.read_structure(args.reference)

    ns = gemmi.NeighborSearch(structure, max_radius=1).populate()

    output = {}

    for c in reference[0]:
        for r in c:
            if (
                r.name not in data
                or gemmi.find_tabulated_residue(r.name).is_amino_acid()
            ):
                continue

            centroid = np.mean([a.pos.tolist() for a in r], axis=0)

            nearest_atom = ns.find_nearest_atom(
                gemmi.Position(centroid[0], centroid[1], centroid[2])
            )
            nearest_residue = structure[0][nearest_atom.chain_idx][
                nearest_atom.residue_idx
            ]
            nearest_chain = structure[0][nearest_atom.chain_idx]

            n = 0
            d_sum = 0

            for reference_atom in r:
                atom = nearest_residue.find_atom(reference_atom.name, "\0")
                if not atom:
                    continue
                if atom.name not in ["C1", "C2", "C3", "C4", "C5", "O5"]:  # ring atoms
                    continue

                d = (reference_atom.pos - atom.pos).length() ** 2
                if d > 2:
                    continue
                d_sum += d
                n += 1

            key = format_residue(c, r)
            alt_key = format_residue(nearest_chain, nearest_residue)
            if n > 0:
                rmsd = np.sqrt(d_sum / n)
                output[key] = {"residue_key": alt_key, "rmsd": round(rmsd, 2)}
            else:
                output[key] = {"residue_key": None, "rmsd": None}

    total_sugars = 0
    modelled = 0
    for k, v in output.items():
        rk = v["residue_key"]
        rrmsd = v["rmsd"]
        if not rk and not rrmsd:
            total_sugars += 1
        else:
            modelled += 1
            total_sugars += 1

    percentage_modelled = 100 * modelled / total_sugars
    print(f"Percentage Modelled {percentage_modelled:.2f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-model", type=str, required=True)
    parser.add_argument("-reference", type=str, required=True)

    args = parser.parse_args()

    main(args)
