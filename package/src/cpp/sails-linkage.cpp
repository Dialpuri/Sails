//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../include/sails-linkage.h"
#include "../include/sails-refine.h"


bool Sails::Model::check_entry_in_database(gemmi::Residue residue) {
    if (linkage_database.find(residue.name) == linkage_database.end()) {
        std::cout << "Residue type: " << residue.name << " is not in Sails' Linkage Database" << std::endl;
        return true;
    }
    return false;
}

Sails::Glycan Sails::Model::extend(Glycan &glycan, int base_seqid, Density &density) {
    const std::vector<Sugar *> terminal_sugars = glycan.get_terminal_sugars(base_seqid);

    for (auto &terminal_sugar: terminal_sugars) {
        auto residue = Utils::get_residue_from_glycosite(terminal_sugar->site, structure);
        // std::cout << "Terminal sugar is " << Utils::format_residue_key(&residue) << std::endl;
        if (check_entry_in_database(residue)) continue;

        for (auto &data: linkage_database[residue.name]) {
            std::optional<SuperpositionResult> opt_result = add_residue(residue, data, density, true);
            if (!opt_result.has_value()) { std::cout << "Could not position another sugar" << std::endl; continue; };

            SuperpositionResult result = opt_result.value();
            result.new_residue.seqid = gemmi::SeqId(terminal_sugar->seqId + 100, 0);

            auto all_residues = &structure.models[terminal_sugar->site.model_idx].chains[terminal_sugar->site.chain_idx]
                    .residues;

            all_residues->insert(all_residues->end(), std::move(result.new_residue));
        }
    }

    Topology topology = {structure, residue_database};
    return topology.find_glycan_topology(glycan.glycosite);

}

void Sails::Model::move_acceptor_atomic_positions(std::vector<gemmi::Atom *> &atoms, double length,
                                                  std::vector<double> &angles, std::vector<double> &torsions) {
    for (int i = 0; i < 3; i++) {
        gemmi::Position atom1 = atoms[i]->pos;
        gemmi::Position atom2 = atoms[i + 1]->pos;
        gemmi::Position atom3 = atoms[i + 2]->pos;
        gemmi::Vec3 new_position = calculate_projected_point(atom1, atom2, atom3, length, angles[i], torsions[i]);
        atoms[i + 3]->pos = gemmi::Position(new_position);
    }
}

gemmi::Transform Sails::Model::superpose_atoms(std::vector<gemmi::Atom *> &atoms,
                                               std::vector<gemmi::Atom> &reference_atoms, double length,
                                               std::vector<double> &angles, std::vector<double> &torsions) {
    move_acceptor_atomic_positions(atoms, length, angles, torsions);

    std::vector<gemmi::Position> reference_positions;
    std::vector<gemmi::Position> new_positions;
    std::for_each(atoms.begin() + 3, atoms.end(), [&](const gemmi::Atom *a) { new_positions.emplace_back(a->pos); });
    std::for_each(reference_atoms.begin(), reference_atoms.end(), [&](const gemmi::Atom &a) {
        reference_positions.emplace_back(a.pos);
    });

    return calculate_superposition(reference_positions, new_positions);
}

void Sails::Model::save_pdb(const std::string &path) {
    std::ofstream of;
    of.open(path);
    write_pdb(structure, of);
    of.close();
}

void Sails::Model::remove_leaving_atom(Sails::LinkageData &data, gemmi::Residue& reference_library_monomer, gemmi::Residue& new_monomer) {
    std::string leaving_atom = "O" + std::to_string(data.acceptor_number); // Ox atom leaves in sugars
    auto leaving_atom_iter = std::remove_if(new_monomer.atoms.begin(), new_monomer.atoms.end(),
                                            [&leaving_atom](const gemmi::Atom &a) { return a.name == leaving_atom; });
    new_monomer.atoms.erase(leaving_atom_iter, new_monomer.atoms.end());

    leaving_atom_iter = std::remove_if(reference_library_monomer.atoms.begin(), reference_library_monomer.atoms.end(),
                                       [&leaving_atom](const gemmi::Atom &a) { return a.name == leaving_atom; });
    reference_library_monomer.atoms.erase(leaving_atom_iter, reference_library_monomer.atoms.end());
}

std::optional<Sails::SuperpositionResult> Sails::Model::add_residue(
    const gemmi::Residue &residue, LinkageData &data, Density &density, bool refine) {
    std::vector<gemmi::Atom *> atoms;

    // find library monomer for acceptor residue
    auto library_monomer = get_monomer(data.acceptor);

    if (!library_monomer.has_value()) { throw std::runtime_error("Could not get monomer"); }
    auto reference_library_monomer = gemmi::Residue(library_monomer.value());
    auto new_monomer = gemmi::Residue(library_monomer.value());

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

    std::vector<double> torsions = data.torsions.get_means_in_order();
    std::vector<double> angles = data.angles.get_in_order();
    double length = data.length;

    gemmi::Transform superpose_result = superpose_atoms(atoms, reference_atoms, length, angles, torsions);
    gemmi::transform_pos_and_adp(new_monomer, superpose_result);

    // remove leaving atom
    remove_leaving_atom(data, reference_library_monomer, new_monomer);

    SuperpositionResult result = {new_monomer, superpose_result, reference_library_monomer};

    if (refine) {
        TorsionAngleRefiner refiner = {atoms, reference_atoms, density, result, length};
        SuperpositionResult final_result = refiner.refine(angles, torsions);
        return final_result;
    }

    float rscc = density.rscc_score(result);
    std::cout << "RSCC = " << rscc << std::endl;
    if (rscc < 0.3) {
        std::cout << "RSCC is less than 0.3" << std::endl;
        return std::nullopt;
    }


    return result;
}


std::optional<gemmi::Residue> Sails::Model::get_monomer(const std::string &monomer) {
    // lookup mononer in cache
    if (monomers.find(monomer) != monomers.end()) {
        return monomers[monomer];
    }

    std::string path = monomer_library_path + "/" + char(std::tolower(monomer.front())) + "/" + monomer + ".cif";
    if (!std::filesystem::exists(path)) {
        std::cout << path << " monomer does not exist" << std::endl;
        return std::nullopt;
    }
    gemmi::Structure structure = gemmi::read_structure_file(path, gemmi::CoorFormat::Detect);
    gemmi::Residue residue = structure.models[0].chains[0].residues[0];

    monomers[monomer] = residue;
    return residue;
}

