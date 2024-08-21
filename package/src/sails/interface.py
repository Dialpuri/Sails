import logging
from typing import Tuple
from pathlib import Path
import gemmi
import sails
import numpy as np


def read_sf_cif(mtz: Path) -> gemmi.Mtz:
    """
    :param mtz: Path to the CIF file containing structure factor data.
    :return: MTZ file containing the converted structure factor data.
    """
    doc = gemmi.cif.read(str(mtz))
    rblocks = gemmi.as_refln_blocks(doc)
    if not rblocks:
        raise RuntimeError("SF CIF supplied has no reflections")
    cif2mtz = gemmi.CifToMtz()
    return cif2mtz.convert_block_to_mtz(rblocks[0])


def get_sails_map(map: gemmi.Ccp4Map | gemmi.FloatGrid | Path | str) -> sails.Grid:
    """
    :param map: The input map, which can be of type gemmi.Ccp4Map, gemmi.FloatGrid, Path or str.
    :return: The Sails map extracted from the input.
    """
    if isinstance(map, gemmi.Ccp4Map):
        sails_grid = extract_gemmi_grid(map.grid)
    elif isinstance(map, gemmi.FloatGrid):
        sails_grid = extract_gemmi_grid(map)
    elif isinstance(map, Path) or isinstance(map, str):
        grid = gemmi.read_ccp4_map(str(map))
        sails_grid = extract_gemmi_grid(grid.grid)
    else:
        raise RuntimeError(
            "Unknown object passed to second argument of n_glycosylate function, allowed types are"
            "gemmi.Mtz, Path, and str"
        )
    return sails_grid


def get_sails_mtz(
    mtz: gemmi.Mtz | Path | str, f: str, sigf: str, fwt: str, phwt: str
) -> sails.MTZ:
    """
    :param mtz: Path to an MTZ file, an instance of gemmi.Mtz, or a string representing the path to an MTZ file.
    :param f: Column name of the F values in the MTZ file.
    :param sigf: Column name of the SIGF values in the MTZ file.
    :param fwt: Column name of the FWT values in the MTZ file.
    :param phwt: Column name of the PHWT values in the MTZ file.
    :return: A gemmi.Mtz object containing only the specified F and SIGF columns.

    This method extracts the specified F and SIGF columns from an MTZ file and returns a new MTZ object
    containing only those columns. It supports three types of input for the 'mtz' parameter:
    1) A string or Path object representing the path to the MTZ file.
    2) An instance of the gemmi.Mtz class.
    3) A gemmi.Mtz object.

    If 'mtz' is an instance of the gemmi.Mtz class, the method directly extracts the specified columns
    using the provided 'f' and 'sigf' column names.

    If 'mtz' is a string or Path object representing the path to an MTZ file, the method first checks
    the file extension to determine the format of the file. If the file has a CIF or ENT extension,
    it is assumed to be a CIF file and the 'read_sf_cif' function is used to read the data. Otherwise,
    the 'gemmi.read_mtz_file' function is used to read the MTZ file. The extracted columns are then
    obtained from the resulting MTZ object using the specified 'f' and 'sigf' column names.

    If 'mtz' is neither an instance of gemmi.Mtz, nor a string/Path object, a RuntimeError is raised.

    The method returns a new gemmi.Mtz object containing only the specified F and SIGF columns.
    """
    if isinstance(mtz, gemmi.Mtz):
        sails_mtz = extract_gemmi_mtz(mtz=mtz, column_names=[f, sigf, fwt, phwt])
    elif isinstance(mtz, Path) or isinstance(mtz, str):
        mtz = Path(mtz)
        if ".cif" in mtz.suffixes or ".ent" in mtz.suffixes:
            m = read_sf_cif(mtz)
        else:
            m = gemmi.read_mtz_file(str(mtz))
        sails_mtz = extract_gemmi_mtz(mtz=m, column_names=[f, sigf, fwt, phwt])
    else:
        raise RuntimeError(
            "Unknown object passed to second argument of n_glycosylate function, allowed types are"
            "gemmi.Mtz, Path, and str"
        )
    return sails_mtz


def get_sails_structure(structure) -> sails.Structure:
    """
    Retrieves the Sails structure from the provided input.

    :param structure: The input structure. It can be either of type gemmi.Structure, or a file path as type Path
                      or str.
    :return: The Sails structure object.
    :raises RuntimeError: If an unknown object is passed as the first argument. Only gemmi.Structure, Path, and str
                          are allowed types.
    """
    if isinstance(structure, gemmi.Structure):
        sails_structure = extract_gemmi_structure(structure=structure)
    elif isinstance(structure, Path) or isinstance(structure, str):
        s = gemmi.read_structure(str(structure))
        sails_structure = extract_gemmi_structure(structure=s)
    else:
        raise RuntimeError(
            "Unknown object passed to first argument of n_glycosylate function, allowed types are"
            "gemmi.Structure, Path, and str"
        )
    return sails_structure


def extract_gemmi_structure(structure: gemmi.Structure) -> sails.Structure:
    os = sails.Structure()
    os.set_cell(
        sails.Cell(
            structure.cell.a,
            structure.cell.b,
            structure.cell.c,
            structure.cell.alpha,
            structure.cell.beta,
            structure.cell.gamma,
        )
    )
    om = sails.Model()
    om.name = structure[0].name
    for chain in structure[0]:
        oc = sails.Chain()
        oc.name = chain.name
        for residue in chain:
            or_ = sails.Residue()
            or_.name = residue.name
            or_.seqid = sails.SeqId(residue.seqid.num, residue.seqid.icode)
            or_.subchain = residue.subchain
            if residue.label_seq:
                or_.set_label_seq(residue.label_seq)

            or_.entity_id = residue.entity_id
            for atom in residue:
                oa = sails.Atom()
                oa.pos = sails.Position(*atom.pos.tolist())
                oa.occ = atom.occ
                oa.b_iso = atom.b_iso
                oa.element = sails.Element(atom.element.name)
                oa.name = atom.name
                oa.set_altloc(atom.altloc)
                or_.atoms.append(oa)
            oc.residues.append(or_)
        om.chains.append(oc)
    os.models.append(om)
    return os


def extract_sails_structure(structure: sails.Structure) -> gemmi.Structure:
    os = gemmi.Structure()
    om = gemmi.Model("1")
    for model in structure.models:
        for chain in model.chains:
            oc = gemmi.Chain(chain.name)
            for residue in chain.residues:
                or_ = gemmi.Residue()
                or_.name = residue.name
                or_.seqid = gemmi.SeqId(residue.seqid.num(), residue.seqid.icode())
                or_.label_seq = residue.get_label_seq()
                or_.entity_id = residue.entity_id
                or_.subchain = residue.subchain
                for atom in residue.atoms:
                    oa = gemmi.Atom()
                    oa.pos = gemmi.Position(atom.pos.x, atom.pos.y, atom.pos.z)
                    oa.occ = atom.occ
                    oa.b_iso = atom.b_iso
                    oa.element = gemmi.Element(atom.element.name)
                    oa.name = atom.name
                    oa.altloc = atom.altloc
                    or_.add_atom(oa)
                oc.add_residue(or_)
            om.add_chain(oc)
        os.add_model(om)

    cell = structure.cell()
    os.cell = gemmi.UnitCell(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)

    return os


def find_alternate_column_labels(mtz: gemmi.Mtz) -> Tuple[str, str]:
    f_type = mtz.columns_with_type("F")
    q_type = mtz.columns_with_type("Q")
    if f_type and q_type:
        best_f_label = f_type[0].label
        best_q_label = q_type[0].label
        if len(q_type) > 1:
            for label in q_type:
                if best_f_label in label.label:
                    best_q_label = label.label
        return best_f_label, best_q_label


def extract_gemmi_mtz(mtz: gemmi.Mtz, column_names=None) -> sails.MTZ:
    if column_names is None:
        column_names = ["FP", "SIGFP", None, None]

    # Find another suitable pair of F,SIGF values if they are not presented
    mtz_labels = mtz.column_labels()
    for name in column_names[:2]:  # only for F and SIGF columns
        if name not in mtz_labels:
            alternate_labels = find_alternate_column_labels(mtz)
            if alternate_labels:
                logging.warning(
                    f"Sails has located two columns of type F and Q ({','.join(alternate_labels)}), check they are "
                    f"correct. Sails will ignore any map coefficient columns."
                )

                column_names = [
                    *alternate_labels,
                    None,
                    None,
                ]  # added None as FWT and PHWT columns
                break
            raise RuntimeError(
                "Could not find suitable column labels in MTZ to continue. Sails requires F and SIGF."
            )

    f, sigf, fwt, phwt = column_names

    fobs = mtz.get_value_sigma(f, sigf)

    fwtobs = {}
    if fwt and phwt:
        f_phi = mtz.get_f_phi(fwt, phwt)
        for hkl, val in zip(f_phi.miller_array, f_phi.value_array):
            fwtobs[tuple(hkl)] = val

    data = []
    for fo in fobs:
        hkl = sails.HKL(*fo.hkl)
        f_sigf = sails.Pair(fo.value.value, fo.value.sigma)
        hkl_key = tuple(fo.hkl)
        if hkl_key in fwtobs:
            fw = fwtobs[hkl_key]
            fwt_phwt = sails.Pair(np.abs(fw), np.rad2deg(np.angle(fw)))
            reflection = sails.Reflection(hkl, f_sigf, fwt_phwt)
        else:
            reflection = sails.Reflection(hkl, f_sigf)

        data.append(reflection)

    cell = sails.Cell(
        mtz.cell.a,
        mtz.cell.b,
        mtz.cell.c,
        mtz.cell.alpha,
        mtz.cell.beta,
        mtz.cell.gamma,
    )
    return sails.MTZ(data, cell, mtz.spacegroup.hm)


def extract_sails_mtz(mtz: sails.MTZ) -> gemmi.Mtz:
    mtz_data = []
    for refln in mtz.reflections:
        data = []
        data.append(refln.hkl.h)
        data.append(refln.hkl.k)
        data.append(refln.hkl.l)
        data.append(refln.f_sigf.value1)
        data.append(refln.f_sigf.value2)
        data.append(refln.fwt_phwt.value1)
        data.append(refln.fwt_phwt.value2)
        data.append(refln.delfwt_phdelwt.value1)
        data.append(refln.delfwt_phdelwt.value2)
        mtz_data.append(data)

    new_mtz = gemmi.Mtz(with_base=True)
    new_mtz.spacegroup = gemmi.SpaceGroup(mtz.spacegroup)
    new_mtz.set_cell_for_all(
        gemmi.UnitCell(
            mtz.cell.a,
            mtz.cell.b,
            mtz.cell.c,
            mtz.cell.alpha,
            mtz.cell.beta,
            mtz.cell.gamma,
        )
    )
    new_mtz.add_dataset("sails")
    new_mtz.add_column("FP", "F")
    new_mtz.add_column("SIGFP", "Q")
    new_mtz.add_column("FWT", "F")
    new_mtz.add_column("PHWT", "P")
    new_mtz.add_column("DELFWT", "F")
    new_mtz.add_column("PHDELWT", "P")
    new_mtz.set_data(np.array(mtz_data))
    new_mtz.ensure_asu()
    new_mtz.update_reso()
    return new_mtz


def extract_gemmi_grid(grid: gemmi.FloatGrid) -> sails.Grid:
    sails_grid = sails.Grid()
    sails_grid.set_data(grid.array)
    sails_grid.nu = grid.nu
    sails_grid.nv = grid.nv
    sails_grid.nw = grid.nw
    sails_grid.set_spacegroup(grid.spacegroup.hm)
    cell = sails.Cell(
        grid.unit_cell.a,
        grid.unit_cell.b,
        grid.unit_cell.c,
        grid.unit_cell.alpha,
        grid.unit_cell.beta,
        grid.unit_cell.gamma,
    )
    sails_grid.set_cell(cell)
    axis_order = sails.AxisOrder(value=grid.axis_order.value)
    sails_grid.set_axis_order(axis_order)
    sails_grid.set_spacing(*grid.spacing)
    return sails_grid
