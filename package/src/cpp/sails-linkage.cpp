//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../include/sails-linkage.h"

#include <gemmi/mmread.hpp>
#include <gemmi/modify.hpp>
#include <gemmi/qcp.hpp>
#include <src/include/sails-topology.h>

Sails::Glycan Sails::Model::extend(Sails::Glycan &glycan, int base_seqid) {
    std::vector<Sugar *> terminal_sugars = glycan.get_terminal_sugars(base_seqid);

    for (const auto &x: terminal_sugars) {
        auto residue = Sails::Utils::get_residue_from_glycosite(x->site, structure);

        if (linkage_database.find(residue.name) == linkage_database.end()) {
            //            std::cout << "Residue type: " << residue.name << " is not in Sails' Linkage Database" << std::endl;
            continue;
        }

        std::cout << residue.name << std::endl;
        for (auto &data: linkage_database[residue.name]) {
            // std::cout << "Attempting " << Sails::Utils::linkage_to_id(data) << std::endl;
            gemmi::Residue new_sugar = calculate_residue_translation(residue, data);
            structure.models[x->site.model_idx].chains[x->site.chain_idx].append_residues({new_sugar});
        }
    }

    Sails::Topology topology = {structure, residue_database};
    if (auto new_glycan = topology.find_glycan_topology(glycan.glycosite); new_glycan.has_value()) return new_glycan.value();
    return {};
}

gemmi::Residue Sails::Model::calculate_residue_translation(const gemmi::Residue &residue, LinkageData &data) {
    std::vector<gemmi::Atom *> atoms;

    // find library monomer for acceptor residue
    auto library_monomer = get_monomer(data.acceptor);

    if (!library_monomer.has_value()) { throw std::runtime_error("Could not get monomer"); }
    auto reference_library_monomer = library_monomer.value();
    std::vector<gemmi::Atom> reference_atoms;

    ResidueData donor_residue = residue_database[data.donor];
    ResidueData acceptor_residue = residue_database[data.acceptor];

    // find donor atoms and add them to atoms list
    auto donor_atoms = donor_residue.donor_map[data.donor_number];
    for (const auto &donor: donor_atoms) {
        auto atom = const_cast<gemmi::Atom *>(residue.find_atom(donor, '*'));
        if (atom == nullptr) { throw std::runtime_error("Could not find an atom"); }
        atoms.emplace_back(atom);
    }

    // find acceptor atoms and add them to atoms list
    auto acceptor_atoms = acceptor_residue.acceptor_map[data.acceptor_number];
    for (const auto &acceptor: acceptor_atoms) {
        auto atom = const_cast<gemmi::Atom *>(library_monomer.value().find_atom(acceptor, '\0'));
        if (atom == nullptr) { throw std::runtime_error("Could not find an atom"); }
        atoms.emplace_back(atom);
        reference_atoms.emplace_back(*atom); // make a copy so that we have a unchanged reference
    }

    if (atoms.size() != 6) { throw std::runtime_error("Unexpected atom count"); }

    TorsionSet linkage_torsions = data.torsions;
    std::vector<double> mean_torsions = data.torsions.get_means_in_order();
    std::vector<double> angles = data.angles.get_in_order();

    for (int i = 0; i < 3; i++) {
        gemmi::Position atom1 = atoms[i]->pos;
        gemmi::Position atom2 = atoms[i + 1]->pos;
        gemmi::Position atom3 = atoms[i + 2]->pos;

        gemmi::Vec3 new_position = calculate_projected_point(atom1, atom2, atom3, 2.0, angles[i], mean_torsions[i]);
        atoms[i + 3]->pos = gemmi::Position(new_position);
    }

    std::vector<gemmi::Position> reference_positions;
    std::vector<gemmi::Position> new_positions;

    std::for_each(atoms.begin() + 3, atoms.end(), [&](const gemmi::Atom *a) { new_positions.emplace_back(a->pos); });

    std::for_each(reference_atoms.begin(), reference_atoms.end(), [&](const gemmi::Atom a) {
        reference_positions.emplace_back(a.pos);
    });

    auto superpose_result = gemmi::superpose_positions(new_positions.data(), reference_positions.data(), 3, nullptr);
    transform_pos_and_adp(reference_library_monomer, superpose_result.transform); // gemmi function
    return reference_library_monomer;
}

std::optional<gemmi::Residue> Sails::Model::get_monomer(const std::string &monomer) {
    std::string path = monomer_library_path + "/" + char(std::tolower(monomer.front())) + "/" + monomer + ".cif";
    if (!std::filesystem::exists(path)) {
        std::cout << path << " monomer does not exist" << std::endl;
        return std::nullopt;
    }
    gemmi::Structure structure = gemmi::read_structure_file(path, gemmi::CoorFormat::Detect);
    return structure.models[0].chains[0].residues[0];
}

