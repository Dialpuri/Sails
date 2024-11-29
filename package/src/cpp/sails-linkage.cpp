//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../include/sails-linkage.h"
#include "../include/sails-refine.h"

// LOGGING FUNCTIONS

void Sails::Model::print_addition_log(const Sails::Sugar *terminal_sugar, Sails::SuperpositionResult &addition) {
    std::cout << "Adding " << Utils::format_residue_key(&addition.new_residue) << " to " <<
            Utils::format_residue_from_site(terminal_sugar->site, structure) << std::endl;
}

void Sails::Model::print_attempted_addition_log(gemmi::Residue *residue_ptr, Sails::LinkageData &data,
                                                const Glycosite *site) {
    std::cout << "\tAttempting to add " << Utils::format_residue_from_site(*site, structure) << "(" << residue_ptr->
            seqid.str() << ")-" << data.
            donor_number << "," << data.acceptor_number << "-" << data.acceptor << "(?)...";
}

void Sails::Model::print_rejection_log() {
    std::cout << " rejected" << std::endl;
}

void Sails::Model::print_successful_log(Sails::Density &density, std::optional<Sails::SuperpositionResult> opt_result) {
    float rscc = density.rscc_score(opt_result.value());
    std::cout << " added " << Utils::format_residue_key(&opt_result.value().new_residue) << " with RSCC = " <<
            rscc << std::endl;
}


// UTILITY FUNCTIONS
std::optional<gemmi::Residue> Sails::Model::get_monomer(const std::string &monomer, bool remove_h) {
    // lookup mononer in cache
    if (monomers.find(monomer) != monomers.end()) {
        return monomers[monomer];
    }

    std::string path = monomer_library_path + "/" + char(std::tolower(monomer.front())) + "/" + monomer + ".cif";

    if (!Utils::file_exists(path)) {
        std::cerr << "File " << path << " does not exist" << std::endl;
        path = special_monomer_path + "/" + monomer + ".cif";
        if (!Utils::file_exists(path)) {
            std::cout << path << " monomer does not exist" << std::endl;
            return std::nullopt;
        }
    }

    gemmi::Structure structure = gemmi::read_structure_file(path, gemmi::CoorFormat::Detect);
    gemmi::Residue residue = structure.models[0].chains[0].residues[0];

    if (remove_h) {
        auto leaving_atom_iter = std::remove_if(residue.atoms.begin(), residue.atoms.end(),
                                                [&](const gemmi::Atom &a) { return a.element.atomic_number() == 1; });
        residue.atoms.erase(leaving_atom_iter, residue.atoms.end());
    }

    monomers[monomer] = residue;
    return residue;
}

std::optional<gemmi::Residue> Sails::Model::get_monomer_only(const std::string &monomer, bool remove_h) {
    std::string monomer_library_path = Utils::get_environment_variable("CLIBD") + "/monomers";
    std::string path = monomer_library_path + "/" + char(std::tolower(monomer.front())) + "/" + monomer + ".cif";
    if (!Utils::file_exists(path)) {
        return std::nullopt;
    }
    gemmi::Structure structure = gemmi::read_structure_file(path, gemmi::CoorFormat::Detect);
    gemmi::Residue residue = structure.models[0].chains[0].residues[0];

    if (remove_h) {
        auto leaving_atom_iter = std::remove_if(residue.atoms.begin(), residue.atoms.end(),
                                                [&](const gemmi::Atom &a) { return a.element.atomic_number() == 1; });
        residue.atoms.erase(leaving_atom_iter, residue.atoms.end());
    }

    return residue;
}


void Sails::Model::save(const std::string &path, std::vector<LinkRecord> &links) {
    std::ofstream os(path);
    gemmi::cif::Document document = make_mmcif_document(*structure);
    gemmi::cif::Block *block = &document.sole_block();
    auto struct_conn = block->find_or_add("_struct_conn", LinkRecord::tags());

    for (LinkRecord &link: links) {
        struct_conn.append_row(link.labels());
    }

    write_cif_to_stream(os, document);
    os.close();
}

void Sails::Model::save(const std::string &path) {
    std::ofstream os(path);
    write_cif_to_stream(os, make_mmcif_document(*structure));
    os.close();
}

bool Sails::Model::residue_in_database(gemmi::Residue *residue) {
    if (linkage_database.find(residue->name) == linkage_database.end()) {
        std::cout << "Residue type: " << residue->name << " is not in Sails' Linkage Database" << std::endl;
        return false;
    }
    return true;
}


// RESIDUE MANIPULATION FUNCTIONS
template<typename T>
void Sails::Model::move_acceptor_atomic_positions(std::vector<T> &atoms, double length,
                                                  std::vector<double> &angles, std::vector<double> &torsions) {
    for (int i = 0; i < 3; i++) {
        gemmi::Position atom1, atom2, atom3;

        if constexpr (std::is_pointer_v<T>) {
            atom1 = atoms[i]->pos;
            atom2 = atoms[i + 1]->pos;
            atom3 = atoms[i + 2]->pos;
        } else {
            atom1 = atoms[i].pos;
            atom2 = atoms[i + 1].pos;
            atom3 = atoms[i + 2].pos;
        }

        gemmi::Vec3 new_position = calculate_projected_point(atom1, atom2, atom3, length, angles[i], torsions[i]);
        if constexpr (std::is_pointer_v<T>) {
            atoms[i + 3]->pos = gemmi::Position(new_position);
        } else {
            atoms[i + 3].pos = gemmi::Position(new_position);
        }
    }
}

gemmi::Transform Sails::Model::superpose_atoms(std::vector<gemmi::Atom *> &atoms,
                                               std::vector<gemmi::Atom> &reference_atoms, double length,
                                               std::vector<double> &angles,
                                               std::vector<double> &torsions) {
    move_acceptor_atomic_positions(atoms, length, angles, torsions);

    std::vector<gemmi::Position> reference_positions;
    std::vector<gemmi::Position> new_positions;
    std::for_each(atoms.begin() + 3, atoms.end(), [&](const gemmi::Atom *a) {
        new_positions.emplace_back(a->pos);
    });
    std::for_each(reference_atoms.begin(), reference_atoms.end(), [&](const gemmi::Atom &a) {
        reference_positions.emplace_back(a.pos);
    });

    return calculate_superposition(reference_positions, new_positions);
}

gemmi::Transform Sails::Model::superpose_atoms(std::vector<gemmi::Atom> &atoms,
                                               std::vector<gemmi::Atom> &reference_atoms, double length,
                                               std::vector<double> &angles,
                                               std::vector<double> &torsions) {
    move_acceptor_atomic_positions(atoms, length, angles, torsions);

    std::vector<gemmi::Position> reference_positions;
    std::vector<gemmi::Position> new_positions;
    std::for_each(atoms.begin() + 3, atoms.end(), [&](const gemmi::Atom a) { new_positions.emplace_back(a.pos); });
    std::for_each(reference_atoms.begin(), reference_atoms.end(), [&](const gemmi::Atom &a) {
        reference_positions.emplace_back(a.pos);
    });

    return calculate_superposition(reference_positions, new_positions);
}


void Sails::Model::remove_leaving_atom(Sails::LinkageData &data, gemmi::Residue &reference_library_monomer,
                                       gemmi::Residue &new_monomer) {
    std::string leaving_atom = "O" + std::to_string(data.acceptor_number); // Ox atom leaves in sugars
    auto leaving_atom_iter = std::remove_if(new_monomer.atoms.begin(), new_monomer.atoms.end(),
                                            [&leaving_atom](const gemmi::Atom &a) { return a.name == leaving_atom; });
    new_monomer.atoms.erase(leaving_atom_iter, new_monomer.atoms.end());

    leaving_atom_iter = std::remove_if(reference_library_monomer.atoms.begin(), reference_library_monomer.atoms.end(),
                                       [&leaving_atom](const gemmi::Atom &a) { return a.name == leaving_atom; });
    reference_library_monomer.atoms.erase(leaving_atom_iter, reference_library_monomer.atoms.end());
}


void Sails::Model::add_sugar_to_structure(const Sugar *terminal_sugar, SuperpositionResult &favoured_addition,
                                          ChainType &chain_type) {
    int chain_idx = terminal_sugar->site.chain_idx;

    if (chain_type == protein) {
        const size_t last_chain_idx = structure->models[terminal_sugar->site.model_idx].chains.size();
        chain_idx = static_cast<int>(last_chain_idx);
        gemmi::Chain chain = gemmi::Chain("");
        chain.name = Utils::get_next_string(
            structure->models[terminal_sugar->site.model_idx].chains[last_chain_idx - 1].name);
        structure->models[terminal_sugar->site.model_idx].chains.emplace_back(chain);
    }

    auto all_residues = &structure->models[terminal_sugar->site.model_idx].chains[chain_idx].residues;
    favoured_addition.new_residue.seqid = gemmi::SeqId(static_cast<int>(all_residues->size()) + 1, '?');
    all_residues->insert(all_residues->end(), std::move(favoured_addition.new_residue));
}

void Sails::Model::rotate_exocyclic_atoms(gemmi::Residue *residue, std::vector<std::string> &atoms, Density &density) {
    gemmi::Atom *a1 = residue->find_atom(atoms[0], '\0');
    gemmi::Atom *a2 = residue->find_atom(atoms[1], '\0');
    gemmi::Atom *a3 = residue->find_atom(atoms[2], '\0');

    if (a1 == nullptr || a2 == nullptr || a3 == nullptr) throw std::runtime_error("Exocylic atom could not be found");
    // create rotation axis vector from position of atom1 to atom2
    gemmi::Vec3 axis = a2->pos - a1->pos;
    axis = axis.normalized(); // normalize the axis vector

    // calculate components of rotation matrix
    gemmi::Position best_pos;
    float best_score = -1e8;

    for (int i = 0; i < 360; i++) {
        double angle_r = clipper::Util::d2rad(i);
        double c = cos(angle_r);
        double s = sin(angle_r);
        double C = 1 - c;
        double x = axis.x;
        double y = axis.y;
        double z = axis.z;

        // create rotation matrix - Rodrigues Rotation Formula
        gemmi::Mat33 rot_matrix(
            x * x * C + c, x * y * C - z * s, x * z * C + y * s,
            y * x * C + z * s, y * y * C + c, y * z * C - x * s,
            z * x * C - y * s, z * y * C + x * s, z * z * C + c
        );

        // rotate around the axis
        auto new_position = gemmi::Position(rot_matrix.multiply(a3->pos - a1->pos) + a1->pos);
        if (const float score = density.score_position(new_position); score > best_score) {
            best_pos = new_position;
            best_score = score;
        }
    }
    a3->pos = best_pos;
}

Sails::Model::ChainType Sails::Model::find_chain_type(std::vector<Sugar *> sugars) {
    if (sugars.empty()) { return protein; }
    gemmi::Chain *chain = Utils::get_chain_ptr_from_glycosite(sugars[0]->site, structure);
    const bool result = std::all_of(chain->residues.begin(), chain->residues.end(), [&](const gemmi::Residue &residue) {
        if (gemmi::find_tabulated_residue(residue.name).is_amino_acid()) return false;
        if (residue_database.find(residue.name) == residue_database.end()) return false;
        return true;
    });
    return result ? non_protein : protein;
}

double Sails::Model::calculate_clash_score(const SuperpositionResult &result) const {
    constexpr double radius = 1;
    gemmi::NeighborSearch ns = gemmi::NeighborSearch(structure->models[0], structure->cell, radius).populate();
    double clash_score = 0;
    for (auto &atom: result.new_residue.atoms) {
        auto nearest_atoms = ns.find_atoms(atom.pos, '\0', 0, radius);
        clash_score += static_cast<double>(nearest_atoms.size());
    }
    return clash_score;
}

std::optional<Sails::SuperpositionResult> Sails::Model::add_residue(
    gemmi::Residue *residue, LinkageData &data, Density &density, bool refine) {
    // find library monomer for acceptor residue
    auto library_monomer = get_monomer(data.acceptor, true);

    if (!library_monomer.has_value()) {
        throw std::runtime_error("Could not get required monomer, "
            "ensure that CCP4 monomer library is sourced. "
            "If you have a local monomer library, ensure that you "
            "have CLIBD set");
    }
    auto reference_library_monomer = gemmi::Residue(library_monomer.value());

    std::vector<gemmi::Atom> reference_atoms;
    ResidueData donor_residue = residue_database[data.donor];
    ResidueData acceptor_residue = residue_database[data.acceptor];

    // find donor atoms and add them to atoms list
    auto donor_atoms = donor_residue.donor_map[data.donor_number];

    if (const std::string x = donor_atoms[2]; x == "O6") {
        // if the donor is an exocyclic atom, rotate the bond
        rotate_exocyclic_atoms(residue, donor_atoms, density);
    }


    std::vector<gemmi::Atom *> atoms;
    for (const auto &donor: donor_atoms) {
        auto atom = residue->find_atom(donor, '*');
        if (atom == nullptr) { throw std::runtime_error("Could not find an atom"); }
        atoms.emplace_back(atom);
    }

    // find acceptor atoms and add them to atoms list
    auto acceptor_atoms = acceptor_residue.acceptor_map[data.acceptor_number];
    for (const auto &acceptor: acceptor_atoms) {
        auto atom = library_monomer.value().find_atom(acceptor, '\0');
        if (atom == nullptr) { throw std::runtime_error("Could not find an atom"); }
        atoms.emplace_back(atom);
        reference_atoms.emplace_back(*atom); // make a copy so that we have a unchanged reference
    }

    if (atoms.size() != 6) { throw std::runtime_error("Unexpected atom count"); }

    double length = data.length;

    SuperpositionResult best_result;
    float best_rscc = INT_MIN;

    for (auto &cluster: data.clusters) {
        std::vector<double> torsions = cluster.torsions.get_means_in_order();
        std::vector<double> torsion_stddev = cluster.torsions.get_stddev_in_order();
        std::vector<double> angles = cluster.angles.get_means_in_order();
        std::vector<double> angles_stddev = cluster.angles.get_stddev_in_order();

        auto new_monomer = gemmi::Residue(library_monomer.value());
        gemmi::Transform superpose_result = superpose_atoms(atoms, reference_atoms, length, angles, torsions);
        gemmi::transform_pos_and_adp(new_monomer, superpose_result);

        // remove leaving atom
        remove_leaving_atom(data, reference_library_monomer, new_monomer);
        new_monomer.seqid = gemmi::SeqId(residue->seqid.num.value + 1, 0);
        reference_library_monomer.seqid = gemmi::SeqId(residue->seqid.num.value + 1, 0);

        SuperpositionResult result = {new_monomer, superpose_result, reference_library_monomer};

        if (refine) {
            TorsionAngleRefiner refiner = {
                atoms, reference_atoms, density, result, length, angles, angles_stddev, torsions, torsion_stddev
            };
            result = refiner.refine();
        }

        // calculate clash score
        double clash_score = calculate_clash_score(result);
        if (clash_score > 1) {
            continue;
        }

        // calculate rscc
        float rscc = density.rscc_score(result);
        if (rscc < 0.2) {
            continue;
        }
        if (rscc > best_rscc) {
            best_rscc = rscc;
            best_result = result;
        }
    }
    if (best_rscc == INT_MIN) return std::nullopt;
    return best_result;
}

std::optional<Sails::SuperpositionResult> Sails::Model::add_residue(gemmi::Residue *residue, LinkageData &data) {
    auto library_monomer = get_monomer(data.acceptor, true);

    if (!library_monomer.has_value()) {
        throw std::runtime_error("Could not get required monomer, "
            "ensure that CCP4 monomer library is sourced. "
            "If you have a local monomer library, ensure that you "
            "have CLIBD set");
    }
    auto reference_library_monomer = gemmi::Residue(library_monomer.value());

    std::vector<gemmi::Atom> reference_atoms;
    ResidueData donor_residue = residue_database[data.donor];
    ResidueData acceptor_residue = residue_database[data.acceptor];

    // find donor atoms and add them to atoms list
    auto donor_atoms = donor_residue.donor_map[data.donor_number];
    std::vector<gemmi::Atom> atoms;
    for (const auto &donor: donor_atoms) {
        gemmi::Atom *atom = residue->find_atom(donor, '*');
        if (atom == nullptr) { throw std::runtime_error("Could not find an atom"); }
        atoms.emplace_back(*atom);
    }

    // find acceptor atoms and add them to atoms list
    auto acceptor_atoms = acceptor_residue.acceptor_map[data.acceptor_number];
    for (const auto &acceptor: acceptor_atoms) {
        auto atom = library_monomer.value().find_atom(acceptor, '\0');
        if (atom == nullptr) { throw std::runtime_error("Could not find an atom"); }
        atoms.emplace_back(*atom);
        reference_atoms.emplace_back(*atom); // make a copy so that we have a unchanged reference
    }

    if (atoms.size() != 6) { throw std::runtime_error("Unexpected atom count"); }

    double length = data.length;

    SuperpositionResult best_result;
    double best_clash = INT_MAX;

    for (auto &cluster: data.clusters) {
        std::vector<double> torsions = cluster.torsions.get_means_in_order();
        std::vector<double> torsion_stddev = cluster.torsions.get_stddev_in_order();
        std::vector<double> angles = cluster.angles.get_means_in_order();
        std::vector<double> angles_stddev = cluster.angles.get_stddev_in_order();

        auto new_monomer = gemmi::Residue(library_monomer.value());
        gemmi::Transform superpose_result = superpose_atoms(atoms, reference_atoms, length, angles, torsions);
        gemmi::transform_pos_and_adp(new_monomer, superpose_result);

        remove_leaving_atom(data, reference_library_monomer, new_monomer);
        new_monomer.seqid = gemmi::SeqId(std::to_string(residue->seqid.num.value + 1));
        reference_library_monomer.seqid = gemmi::SeqId(std::to_string(residue->seqid.num.value + 1));

        SuperpositionResult result = {new_monomer, superpose_result, reference_library_monomer};

        // calculate clash score
        double clash_score = calculate_clash_score(result);
        if (clash_score < best_clash) {
            best_clash = clash_score;
            best_result = std::move(result);
        }
    }

    if (best_clash == INT_MAX) return std::nullopt;
    return best_result;
}

void Sails::Model::create_pseudo_glycan(PseudoGlycan &pseudo_glycan) {
    // gemmi::Residue *root_residue = Utils::get_residue_ptr_from_glycosite(pseudo_glycan.root_glycosite, structure);
    // auto real_residue = get_monomer(root_residue->name, true);
    // if (!real_residue.has_value()) { throw std::runtime_error("Could not get real residue"); }
    // root_residue->atoms = real_residue.value().atoms;
    typedef std::pair<Sugar *, Sugar *> SugarPair;

    std::stack<SugarPair> sugar_stack;
    std::unordered_set<Sugar *> visited_sugars;

    sugar_stack.emplace(nullptr, pseudo_glycan.root_sugar);

    while (!sugar_stack.empty()) {
        auto [parent, node] = sugar_stack.top();
        sugar_stack.pop();

        if (visited_sugars.find(node) == visited_sugars.end()) {
            if (parent != nullptr) {
                gemmi::Residue *parent_residue = Utils::get_residue_ptr_from_glycosite(parent->site, structure);
                gemmi::Residue *current_residue = Utils::get_residue_ptr_from_glycosite(node->site, structure);

                const Linkage *linkage = parent->find_linkage(parent, node);

                std::cout << Utils::format_residue_from_site(parent->site, structure) << std::endl;

                if (linkage_database.find(parent_residue->name) == linkage_database.end()) {
                    throw std::runtime_error("Could not find entry in linkage database");
                }
                std::vector<LinkageData> available_linkages = linkage_database[parent_residue->name];

                auto found_linkage = std::find_if(available_linkages.begin(), available_linkages.end(),
                                                  [&](const LinkageData &available_linkage) {
                                                      return available_linkage.donor_number == linkage->donor_number &&
                                                             available_linkage.acceptor == current_residue->name;
                                                  });

                if (found_linkage != available_linkages.end()) {
                    std::optional<SuperpositionResult> result = add_residue(parent_residue, *found_linkage);
                    if (!result.has_value()) { throw std::runtime_error("Could not add residue"); }

                    gemmi::Chain *parent_chain = Utils::get_chain_ptr_from_glycosite(node->site, structure);
                    parent_chain->residues[node->site.residue_idx].atoms = result.value().new_residue.atoms;
                } else {
                    std::cout << "Couldn't find a linkage between " << parent_residue->name << " and " <<
                            current_residue->name << std::endl;
                }
            }

            visited_sugars.insert(node);
        }

        for (Sugar *child: pseudo_glycan.adjacency_list[node]) {
            if (visited_sugars.find(child) == visited_sugars.end()) {
                sugar_stack.emplace(node, child);
            }
        }
    }
}

gemmi::Residue Sails::Model::replace_residue(gemmi::Residue *target_residue,
                                             const std::string &replacement_residue_name) {
    const std::vector<std::string> ring_atoms = {"C1", "C2", "C3", "C4", "C5", "O5"};

    auto replacement_residue_result = get_monomer_only(replacement_residue_name, true);
    if (!replacement_residue_result.has_value()) {
        throw std::runtime_error("Could not find replacement residue of specified name - " + replacement_residue_name);
    }

    gemmi::Residue replacement_residue = replacement_residue_result.value();

    std::vector<gemmi::Position> target_positions;
    std::vector<gemmi::Position> source_positions;

    for (auto &ring_atom: ring_atoms) {
        gemmi::Atom *found_target_atom = target_residue->find_atom(ring_atom, '*');
        if (found_target_atom == nullptr) { throw std::runtime_error("Could not find atom in target - " + ring_atom); }

        gemmi::Atom *found_source_atom = replacement_residue.find_atom(ring_atom, '*');
        if (found_source_atom == nullptr) { throw std::runtime_error("Could not find atom in source"); }

        target_positions.emplace_back(found_target_atom->pos);
        source_positions.emplace_back(found_source_atom->pos);
    }

    const gemmi::Transform transform = calculate_superposition(source_positions, target_positions);
    gemmi::transform_pos_and_adp(replacement_residue, transform);
    return replacement_residue;
}


std::vector<Sails::SuperpositionResult> Sails::Model::find_favoured_additions(
    const Sugar *&terminal_sugar, PossibleAdditions &possible_additions) {
    std::vector<SuperpositionResult> favoured_additions;
    for (const auto &[donor, possible_additions]: possible_additions) {
        for (const auto &addition: possible_additions) {
            std::vector<int> preferred_depths = residue_database[addition.new_residue.name].preferred_depths;
            if (std::find(preferred_depths.begin(), preferred_depths.end(), terminal_sugar->depth + 1) !=
                preferred_depths.end()) {
                favoured_additions.emplace_back(addition);
            }
        }
    }
    return favoured_additions;
}


void Sails::Model::extend_if_possible(Density &density, bool debug, ChainType &chain_type,
                                      const Sugar *terminal_sugar) {
    gemmi::Residue *residue_ptr = Utils::get_residue_ptr_from_glycosite(terminal_sugar->site, structure);

    // std::cout << "Terminal sugar is " << Utils::format_site_key(terminal_sugar->site) << " with depth = " << terminal_sugar->depth<< std::endl;
    if (!residue_in_database(residue_ptr)) return;

    // calculate the sugars that can be linked to this terminal sugar
    PossibleAdditions possible_additions;
    for (LinkageData &data: linkage_database[residue_ptr->name]) {
        if (debug) print_attempted_addition_log(residue_ptr, data, &terminal_sugar->site);

        std::optional<SuperpositionResult> opt_result = add_residue(residue_ptr, data, density, true);

        if (!opt_result.has_value()) {
            if (debug) print_rejection_log();
            continue;
        }

        possible_additions[data.donor_number].emplace_back(opt_result.value());

        if (debug) print_successful_log(density, opt_result);
    }

    if (possible_additions.empty()) { return; }

    // find the most likely one for each donor atom  based on preferred depth
    std::vector<SuperpositionResult> favoured_additions = find_favoured_additions(terminal_sugar, possible_additions);

    // add the favoured addition to the structure member
    for (SuperpositionResult &addition: favoured_additions) {
        if (debug) print_addition_log(terminal_sugar, addition);
        add_sugar_to_structure(terminal_sugar, addition, chain_type);
        chain_type = non_protein;
    }
}

Sails::Glycan Sails::Model::extend(Glycan &glycan, Glycosite &base_glycosite, Density &density, bool debug) {
    const std::vector<Sugar *> terminal_sugars = glycan.get_terminal_sugars(base_glycosite);
    ChainType chain_type = find_chain_type(terminal_sugars);

    for (const auto &terminal_sugar: terminal_sugars) {
        extend_if_possible(density, debug, chain_type, terminal_sugar);
    }

    Topology topology = {structure, residue_database};
    return topology.find_glycan_topology(glycan.glycosite);
}
