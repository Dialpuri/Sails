from typing import Tuple

import numpy as np
import sails
import time
import gemmi
from pathlib import Path


def extract_gemmi_structure(structure: gemmi.Structure):
    os = sails.Structure()
    om = sails.Model()
    for chain in structure[0]:
        oc = sails.Chain()
        oc.name = chain.name
        for residue in chain:
            or_ = sails.Residue()
            or_.name = residue.name
            or_.seqid = sails.SeqId(residue.seqid.num, residue.seqid.icode)
            for atom in residue:
                oa = sails.Atom()
                oa.pos = sails.Position(*atom.pos.tolist())
                oa.occ = atom.occ
                oa.b_iso = atom.b_iso
                oa.element = sails.Element(atom.element.name)
                oa.name = atom.name
                or_.atoms.append(oa)
            oc.residues.append(or_)
        om.chains.append(oc)
    os.models.append(om)
    return os


def extract_sails_structure(structure: sails.Structure) -> gemmi.Structure:
    os = gemmi.Structure()
    om = gemmi.Model("")
    for model in structure.models:
        for chain in model.chains:
            oc = gemmi.Chain(chain.name)
            for residue in chain.residues:
                or_ = gemmi.Residue()
                or_.name = residue.name
                or_.seqid = gemmi.SeqId(residue.seqid.num(), residue.seqid.icode())
                for atom in residue.atoms:
                    oa = gemmi.Atom()
                    oa.pos = gemmi.Position(atom.pos.x, atom.pos.y, atom.pos.z)
                    oa.occ = atom.occ
                    oa.b_iso = atom.b_iso
                    oa.element = gemmi.Element(atom.element.name)
                    oa.name = atom.name
                    or_.add_atom(oa)
                oc.add_residue(or_)
            om.add_chain(oc)
        os.add_model(om)

    cell = structure.cell()
    os.cell = gemmi.UnitCell(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)

    return os


def extract_gemmi_mtz(mtz: gemmi.Mtz, column_names=None) -> sails.MTZ:
    if column_names is None:
        column_names = ["FP", "SIGFP"]

    fobs = mtz.get_value_sigma(*column_names)

    data = []
    for entry in fobs:
        f_sigf = sails.Pair(entry.value.value, entry.value.sigma)
        hkl = sails.HKL(*entry.hkl)
        reflection = sails.Reflection(hkl, f_sigf)
        data.append(reflection)

    cell = sails.Cell(mtz.cell.a, mtz.cell.b, mtz.cell.c, mtz.cell.alpha, mtz.cell.beta, mtz.cell.gamma)
    return sails.MTZ(data, cell, mtz.spacegroup_name)


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
        gemmi.UnitCell(mtz.cell.a, mtz.cell.b, mtz.cell.c, mtz.cell.alpha, mtz.cell.beta, mtz.cell.gamma))
    new_mtz.add_dataset('sails')
    new_mtz.add_column('FP', 'F')
    new_mtz.add_column('SIGFP', 'Q')
    new_mtz.add_column('FWT', 'F')
    new_mtz.add_column('PHWT', 'P')
    new_mtz.add_column('DELFWT', 'F')
    new_mtz.add_column('PHDELWT', 'P')
    new_mtz.set_data(np.array(mtz_data))
    new_mtz.ensure_asu()
    new_mtz.update_reso()
    return new_mtz


def n_glycosylate(structure: gemmi.Structure | Path | str, mtz: gemmi.Mtz | Path | str) -> Tuple[
    gemmi.Structure, gemmi.Mtz]:
    sails_structure = None
    sails_mtz = None

    if isinstance(structure, gemmi.Structure):
        sails_structure = extract_gemmi_structure(structure=structure)
    elif isinstance(structure, Path) or isinstance(structure, str):
        s = gemmi.read_structure(str(structure))
        sails_structure = extract_gemmi_structure(structure=s)
    else:
        raise RuntimeError("Unknown object passed to first argument of n_glycosylate function, allowed types are"
                           "gemmi.Structure, Path, and str")

    if isinstance(mtz, gemmi.Mtz):
        sails_mtz = extract_gemmi_mtz(mtz=mtz)
    elif isinstance(mtz, Path) or isinstance(mtz, str):
        m = gemmi.read_mtz_file(str(mtz))
        sails_mtz = extract_gemmi_mtz(mtz=m)
    else:
        raise RuntimeError("Unknown object passed to second argument of n_glycosylate function, allowed types are"
                           "gemmi.Mtz, Path, and str")

    result = sails.n_glycosylate_from_objects(sails_structure, sails_mtz, 1)

    return extract_sails_structure(result.structure), extract_sails_mtz(result.mtz)


def run():
    print("Running Sails")
    t0 = time.time()

    s, m = n_glycosylate("package/models/5fji/5fji.cif", "package/models/5fji/5fji.mtz", )

    s.make_mmcif_block().write_file("sails-5fji.cif")
    m.write_to_file("sails-5fji.mtz")

    t1 = time.time()
    print(f"Sails - Time Taken = {(t1 - t0)} seconds")
