import gemmi
import argparse

def main(args):
    s = gemmi.read_structure(args.pdbin)

    os = gemmi.Structure()
    om = gemmi.Model(s[0].name)
    for c in s[0]:
        oc = gemmi.Chain(c.name)
        for r in c:
            k = gemmi.find_tabulated_residue(r.name)
            if args.rm_glycans:
                if not k.is_amino_acid() and not k.is_nucleic_acid():
                    del r
                    continue
            if args.rm_waters:
                if k.is_water():
                    del r
                    continue
            oc.add_residue(r)
        om.add_chain(oc)
    os.add_model(om)

    if args.pdbin.endswith('.cif'):
        os.make_mmcif_document().write_file(args.pdbout)
    else:
        os.write_pdb(args.pdbout)



if __name__ == "__main__":
    # python deglycosylate_models.py -pdbin ../data/5fji.pdb -pdbout ../data/5fji_dg.pdb -rm-glycans -rm-waters
    parser = argparse.ArgumentParser(
        prog="Deglycosylate"
    )
    parser.add_argument(
        "-pdbin",
        required=True
    )

    parser.add_argument(
        "-pdbout",
        required=True
    )
    parser.add_argument(
        "-rm-glycans",
        action="store_true"
    )
    parser.add_argument(
        "-rm-waters",
        action="store_true"
    )

    args = parser.parse_args()

    main(args)