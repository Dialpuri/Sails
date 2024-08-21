import subprocess
import argparse
import urllib.request

from pathlib import Path


def main(args):
    pdb: str = args.pdb.upper()
    sf_url = f"https://files.rcsb.org/download/{pdb}-sf.cif"
    cif_url = f"https://files.rcsb.org/download/{pdb}.cif"

    pdb_lower = Path(pdb.lower())
    sf_path = pdb_lower / f"{pdb}-sf.cif"
    cif_path = pdb_lower / f"{pdb}.cif"
    mtz_path = pdb_lower / f"{pdb}.mtz"
    deglyco_cif_path = pdb_lower / f"{pdb}_deglycosylated.cif"

    pdb_lower.mkdir(parents=True, exist_ok=True)

    if not sf_path.exists():
        urllib.request.urlretrieve(sf_url, sf_path)

    if not cif_path.exists():
        urllib.request.urlretrieve(cif_url, cif_path)

    if not mtz_path.exists():
        subprocess.run(["gemmi", "cif2mtz", str(sf_path), str(mtz_path)])

    if not deglyco_cif_path.exists():
        subprocess.run(
            [
                "python",
                "/Users/dialpuri/Development/sails/package/scripts/deglycosylate_models.py",
                "-pdbin",
                str(cif_path),
                "-pdbout",
                str(deglyco_cif_path),
                "-rm-glycan",
            ]
        )

    output = pdb_lower / "output"
    output.mkdir(parents=True, exist_ok=True)

    pdb_out = output / f"{pdb.lower()}_sails.cif"
    mtz_out = output / f"{pdb.lower()}_sails.mtz"

    rerun = True
    if not pdb_out.exists() or not mtz_out.exists() or rerun:
        subprocess.run(
            [
                "sails",
                "-pdbin",
                str(deglyco_cif_path),
                "-mtzin",
                str(mtz_path),
                "-pdbout",
                str(pdb_out),
                "-mtzout",
                str(mtz_out),
                "-cycles",
                "6",
            ]
        )

    subprocess.run(
        [
            "python",
            "/Users/dialpuri/Development/sails/package/scripts/compare_structures.py",
            "-model",
            pdb_out,
            "-reference",
            cif_path,
        ]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", type=str, required=True)

    main(parser.parse_args())
