"""
Refine Carbohydrates with CHAPI - Jordan Dialpuri
"""

import argparse
import chapi
import gemmi


def check_carbohydrate(name: str) -> bool:
    """Check if a residue is a carbohydrate using gemmi tables"""
    kind = gemmi.find_tabulated_residue(name)
    return kind.kind == gemmi.ResidueKind.PYR


def refine(args: argparse.Namespace):
    """Refine carbohydrates in a PDB file using CHAPI"""
    mc = chapi.molecules_container_t(False)
    mc.set_make_backups(False)
    mc.set_use_gemmi(False)
    mc.set_refinement_is_verbose(False)
    imol = mc.read_coordinates(args.pdbin)
    if args.mtzin:
        imol_mtz = mc.read_mtz(args.mtzin, "FWT", "PHWT", "W", False, False)
        mc.set_imol_refinement_map(imol_mtz)
    if args.mapin:
        imol_map = mc.read_ccp4_map(args.mapin, False)
        mc.set_imol_refinement_map(imol_map)

    non_standard_codes = mc.non_standard_residue_types_in_model(imol)
    for non_standard_code in non_standard_codes:
        mc.get_monomer_from_dictionary(non_standard_code, imol, False)
        _ = mc.get_monomer(non_standard_code)  # needed for some older chapies.

    structure = gemmi.read_structure(args.pdbin)
    carbohydrate_chains = {}
    for chain in structure[0]:
        carbohydrate_chain = all(
            map(lambda residue: check_carbohydrate(residue.name), chain)
        )
        if not carbohydrate_chain:
            continue

        start_residue = chain[0].seqid.num
        end_residue = chain[-1].seqid.num
        carbohydrate_chains[chain] = (start_residue, end_residue)

    for chain_id, residue_range in carbohydrate_chains.items():
        residue_id_start, residue_id_end = residue_range
        success = mc.refine_residue_range(
            imol, chain_id.name, residue_id_start, residue_id_end, 10000
        )
        if not success:
            print(
                f"Failed to refine Chain {chain_id.name} Residues {residue_id_start}-{residue_id_end}"
            )

    _ = mc.write_coordinates(imol, args.pdbout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Refine carbohydrates")
    parser.add_argument("-pdbin", type=str, required=True)
    parser.add_argument("-mtzin", type=str, required=False)
    parser.add_argument("-mapin", type=str, required=False)
    parser.add_argument("-pdbout", type=str, required=True)
    args = parser.parse_args()

    if not args.mapin and not args.mtzin:
        raise RuntimeError("Must specify either an MTZ or map file")

    refine(args)
