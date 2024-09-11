import argparse
import logging
import os
from sails import interface, get_snfg, get_all_snfgs
import importlib.resources
from pathlib import Path


def main(args):
    if args.all and args.chain:
        logging.warning("--all and -chain were both specified, ignoring chain")

    if args.all and args.seqid:
        logging.warning("--all and -seqid were both specified, ignoring seqid")

    if args.all:
        create_all_snfgs(args)
        return

    if args.chain and args.seqid:
        create_single_snfg(args)
        return

    logging.warning("No valid parameters were specified, no SNFGs have been generated")


def create_all_snfgs(args):
    """
    :param args: A namespace object containing command line arguments
    :return: None


    This method creates all SNFG (Symbol Nomenclature for Glycans) diagrams for the given PDB file.
    It retrieves the Sails structure using the input PDB file, and then reads the required data files from the package.

    If the supplied output path is a file, a directory will be created with the name "sails-snfgs" in the parent directory
    and the SNFG diagrams will be saved in that directory. If the output directory already exists and is empty,
    a fatal error will be logged. Otherwise, the output directory will be used as-is.

    The method uses the get_all_snfgs() function to generate SNFG diagrams for each unique SNFG name in the Sails
    structure. The generated diagrams are saved as SVG files in the output directory.

    Example usage:
        args = Namespace(pdbin='input.cif', svgout='output', "all": True)
        create_all_snfgs(args)
    """
    sails_structure = interface.get_sails_structure(args.model)
    resource = importlib.resources.files("sails").joinpath("data")

    supplied_output_path = Path(args.snfgout)
    output_dir = supplied_output_path

    if not supplied_output_path.is_dir():
        logging.warning("Supplied path is not a directory, one will be created")
        output_dir = supplied_output_path.parent / "sails-snfgs"
        if output_dir.exists() and not next(output_dir.iterdir(), None):
            logging.fatal(
                "Tried to create a directory called snfgs, but one already exists. Please supply an output directory."
            )
            return

        os.makedirs(output_dir, exist_ok=False)

    snfgs = get_all_snfgs(sails_structure, str(resource))
    for k, v in snfgs.items():
        path = output_dir / f"{k}.svg"
        with open(path, "w") as f:
            f.write(v)


def create_single_snfg(args):
    """
    :param args: A namespace object containing command line arguments
    :return: None

    This method takes a dictionary containing command line arguments and uses them to create a specified snfg (sails structure).
    It first retrieves the sails structure from the provided pdbin file using the `get_sails_structure` method from the `interface` module.

    Next, it checks the supplied output path to ensure it does not already exist.
    If it exists, it logs a fatal error and exits the method. If the supplied output path does not have any suffixes
    (i.e., it is a directory), it creates a file inside the directory using the provided chain and seqid. If the output
    file already exists, it logs a fatal error and exits the method.

    Finally, it calls the `get_snfg` method with the chain, seqid, sails structure, and resource path to obtain the
    snfg. It then opens the output file in write mode and writes the snfg content to the file.
    """
    sails_structure = interface.get_sails_structure(args.model)
    resource = importlib.resources.files("sails").joinpath("data")

    supplied_output_path = Path(args.snfgout)
    output_dir = supplied_output_path

    if supplied_output_path.exists() and not args.overwrite:
        logging.fatal(f"Output path already exists {supplied_output_path}. Exiting.")
        return

    if not supplied_output_path.suffixes:
        logging.warning("Supplied path is a directory, file will be created inside")
        supplied_output_path.mkdir(exist_ok=True)
        output_dir = supplied_output_path / f"{args.chain}-{args.seqid}.svg"
        if output_dir.exists():
            logging.fatal(f"Output path already exists {output_dir}. Exiting.")
            return

    snfg = get_snfg(args.chain, args.seqid, sails_structure, str(resource))
    with open(output_dir, "w") as f:
        f.write(snfg)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-model", type=str, required=True)
    parser.add_argument("-snfgout", type=str, required=True)
    parser.add_argument("-chain", type=str, required=False)
    parser.add_argument("-seqid", type=int, required=False)
    parser.add_argument("--all", action=argparse.BooleanOptionalAction, required=False)
    parser.add_argument(
        "--overwrite", action=argparse.BooleanOptionalAction, required=False
    )

    return parser.parse_args()


def run():
    args = parse_args()
    main(args)
