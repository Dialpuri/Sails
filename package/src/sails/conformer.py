import importlib
from pathlib import Path
from sails import get_conformations
import argparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--donor_residue", help="Name of donor residue e.g. NAG", required=True
    )
    parser.add_argument(
        "--acceptor_residue", help="Name of donor residue e.g. NAG", required=True
    )
    parser.add_argument(
        "--donor_number", help="Number of donor atom e.g. 4", required=True
    )
    parser.add_argument(
        "--acceptor_number", help="Number of donor atom e.g. 1", required=True
    )
    parser.add_argument(
        "--step", help="Step size to use for torsion angles, default=6", default=6
    )
    parser.add_argument(
        "--lower_bound",
        help="Upper bound for torsion angle, default=-180",
        default=-180,
    )
    parser.add_argument(
        "--upper_bound", help="Upper bound for torsion angle, default=180", default=180
    )
    parser.add_argument("--output", help="Path to output directory")

    return parser.parse_args()


def run():
    args = parse_args()
    resource = str(importlib.resources.files("sails").joinpath("data"))

    if not args.output:
        output = Path(
            f"{args.donor_residue}-{args.donor_number},{args.acceptor_number}-{args.acceptor_residue}-conformers"
        )
    else:
        output = Path(args.output)

    output.mkdir(exist_ok=True, parents=True)

    get_conformations(
        args.donor_residue,
        args.acceptor_residue,
        args.donor_number,
        args.acceptor_number,
        resource,
        args.step,
        args.lower_bound,
        args.upper_bound,
        str(output),
    )
