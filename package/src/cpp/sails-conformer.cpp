//
// Created by Jordan Dialpuri on 22/09/2024.
//

#include "../include/sails-conformer.h"
#include "src/include/sails-linkage.h"

void Sails::Conformer::remove_leaving_atoms(Sails::LinkageData &data, gemmi::Residue &reference_library_monomer,
                                       gemmi::Residue &new_monomer) {
    std::string leaving_hydrogen = "HO" + std::to_string(data.donor_number); // HOx atom leaves in sugars
    auto leaving_hydrogen_iter = std::remove_if(reference_library_monomer.atoms.begin(), reference_library_monomer.atoms.end(),
                                            [&leaving_hydrogen](const gemmi::Atom &a) { return a.name == leaving_hydrogen; });
    if (leaving_hydrogen_iter != reference_library_monomer.atoms.end()) {
        reference_library_monomer.atoms.erase(leaving_hydrogen_iter, reference_library_monomer.atoms.end());
    }

    std::string leaving_atom = "O" + std::to_string(data.donor_number); // Ox atom leaves in sugars
    auto leaving_atom_iter = std::remove_if(reference_library_monomer.atoms.begin(), reference_library_monomer.atoms.end(),
                                            [&leaving_atom](const gemmi::Atom &a) { return a.name == leaving_atom; });
    reference_library_monomer.atoms.erase(leaving_atom_iter, reference_library_monomer.atoms.end());


    std::string leaving_hydrogen2 = "HO" + std::to_string(data.acceptor_number); // HOx atom leaves in sugars
    auto leaving_hydrogen2_iter = std::remove_if(new_monomer.atoms.begin(), new_monomer.atoms.end(),
                                            [&leaving_hydrogen2](const gemmi::Atom &a) { return a.name == leaving_hydrogen2; });
    if (leaving_atom_iter != new_monomer.atoms.end()) {
        new_monomer.atoms.erase(leaving_hydrogen2_iter, new_monomer.atoms.end());
    }

    // leaving_atom_iter = std::remove_if(reference_library_monomer.atoms.begin(), reference_library_monomer.atoms.end(),
    //                                    [&leaving_atom](const gemmi::Atom &a) { return a.name == leaving_atom; });
    // reference_library_monomer.atoms.erase(leaving_atom_iter, reference_library_monomer.atoms.end());
}

void Sails::Conformer::superimpose(const std::string &residue_1_name, const std::string &residue_2_name,
                                   int donor_number, int acceptor_number, double phi, double psi) {
    std::optional<gemmi::Residue> residue_1 = get_monomer(residue_1_name, false);
    std::optional<gemmi::Residue> residue_2 = get_monomer(residue_2_name, false);
    gemmi::Residue acceptor = residue_2.value();

    std::vector<LinkageData> data = linkage_database[residue_1_name];
    auto it = std::find_if(data.begin(), data.end(), [&](const LinkageData &entry) {
        return entry.donor_number == donor_number && entry.acceptor_number == acceptor_number && entry.acceptor ==
               residue_2_name;
    });

    if (it == data.end()) {
        throw std::runtime_error("Could not find linkage in database");
    }

    ResidueData donor_residue = residue_database[residue_1_name];
    ResidueData acceptor_residue = residue_database[residue_2_name];

    auto donor_atoms = donor_residue.donor_map[donor_number];
    std::vector<gemmi::Atom> reference_atoms;

    std::vector<gemmi::Atom *> atoms;
    for (const auto &donor: donor_atoms) {
        auto atom = residue_1->find_atom(donor, '*');
        if (atom == nullptr) { throw std::runtime_error("Could not find an atom"); }
        atoms.emplace_back(atom);
    }

    // find acceptor atoms and add them to atoms list
    auto acceptor_atoms = acceptor_residue.acceptor_map[acceptor_number];
    for (const auto &acceptor: acceptor_atoms) {
        auto atom = residue_2.value().find_atom(acceptor, '\0');
        if (atom == nullptr) { throw std::runtime_error("Could not find an atom"); }
        atoms.emplace_back(atom);
        reference_atoms.emplace_back(*atom); // make a copy so that we have a unchanged reference
    }

    if (atoms.size() != 6) { throw std::runtime_error("Unexpected atom count"); }

    double length = 1.4;

    std::vector<double> angles = it->clusters[0].angles.get_means_in_order();
    std::vector<double> mean_torsions = it->clusters[0].torsions.get_means_in_order();

    std::vector<double> torsions = {phi, psi, mean_torsions[2]};
    gemmi::Transform superpose_result = superpose_atoms(atoms, reference_atoms, length, angles, torsions);
    gemmi::transform_pos_and_adp(acceptor, superpose_result);
    remove_leaving_atoms(*it, residue_1.value(), acceptor);

    auto x = {residue_1.value(), acceptor};
    std::string path = "../debug/conformers/phi=" + std::to_string(int(phi)) + "_psi=" + std::to_string(int(psi)) + ".pdb";
    Utils::save_residues_to_file(x, path);
}

void Sails::Conformer::get_all_conformers(const std::string &residue_1_name, const std::string &residue_2_name,
                                          int donor_number, int acceptor_number, int step) {
    for (int i = -180; i <= 180; i += step) {
        for (int j = -180; j <= 180; j += step) {
            superimpose(residue_1_name, residue_2_name, donor_number, acceptor_number, (double)i, (double)j);
        }
    }
}
