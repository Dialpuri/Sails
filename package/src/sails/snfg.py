import argparse
import logging
import os

from sails import interface, get_snfg, get_all_snfgs
import importlib.resources
from pathlib import Path
def main(args):

    if args.all:
        create_all_snfgs(args)
        return

    if args.chain and args.seqid:
        create_single_snfg(args)

def create_all_snfgs(args):
    sails_structure = interface.get_sails_structure(args.pdbin)
    resource = importlib.resources.files('sails').joinpath("data")

    supplied_output_path = Path(args.svgout)
    output_dir = supplied_output_path

    if not supplied_output_path.is_dir():
        logging.warning("Supplied path is not a directory, one will be created")
        output_dir = supplied_output_path.parent / "sails-snfgs"
        if output_dir.exists() and not next(output_dir.iterdir(), None):
            logging.fatal(
                "Tried to create a directory called snfgs, but one already exists. Please supply an output directory.")
            return

        os.makedirs(output_dir, exist_ok=False)

    snfgs = get_all_snfgs(sails_structure, str(resource))
    for k, v in snfgs.items():
        path = output_dir / f"{k}.svg"
        with open(path, "w") as f:
            f.write(v)

def create_single_snfg(args):
    sails_structure = interface.get_sails_structure(args.pdbin)
    resource = importlib.resources.files('sails').joinpath("data")

    supplied_output_path = Path(args.svgout)
    output_dir = supplied_output_path

    if supplied_output_path.exists():
        logging.fatal(
            f"Output path already exists {supplied_output_path}. Exiting.")
        return

    if not supplied_output_path.suffixes:
        logging.warning("Supplied path is a directory, file will be created inside")
        supplied_output_path.mkdir(exist_ok=True)
        output_dir = supplied_output_path / f"{args.chain}-{args.seqid}.svg"
        if output_dir.exists():
            logging.fatal(
                f"Output path already exists {output_dir}. Exiting.")
            return

    snfg = get_snfg(args.chain, args.seqid, sails_structure, str(resource))
    with open(output_dir, "w") as f:
        f.write(snfg)

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-pdbin", type=str, required=True)
    parser.add_argument("-svgout", type=str,  required=True)
    parser.add_argument("-chain", type=str, required=False)
    parser.add_argument("-seqid", type=int, required=False)
    parser.add_argument("--all", action=argparse.BooleanOptionalAction, required=False)

    return parser.parse_args()


def run():
    args = parse_args()
    main(args)