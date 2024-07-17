//
// Created by Jordan Dialpuri on 06/07/2024.
//

#include "../include/sails-topology.h"


void Sails::Topology::initialise_neighbour_search(gemmi::Structure* structure) {
    constexpr double search_radius = 2;
    gemmi::NeighborSearch neighbor_search = {structure->models[0], structure->cell, search_radius};
    neighbor_search.populate(false);
    m_neighbor_search = neighbor_search;
}

Sails::Topology::Topology(gemmi::Structure* structure, ResidueDatabase& database) {
    m_structure = structure;
    m_database = database;

    initialise_neighbour_search(structure);
}

void Sails::Topology::find_residue_near_donor(Glycosite &glycosite, Glycan &glycan, std::queue<Glycosite> &queue) {

    double search_radius = 2.5;
    gemmi::Residue residue = Utils::get_residue_from_glycosite(glycosite, m_structure);

    if (m_database.find(residue.name) == m_database.end()) { throw std::runtime_error("Glycosite is not in database"); }
    auto database_entry = m_database[residue.name];

    for (const auto &donor: database_entry.donors) {
        // get donor atoms with that name, could return > 1 with altconfs
        gemmi::AtomGroup donor_atoms = residue.get(donor.atom3);
        for (const auto &donor_atom: donor_atoms) {

            glycan.add_sugar(donor_atom.name, residue.seqid.num.value, glycosite);

            auto near_atoms = m_neighbor_search.find_atoms(donor_atom.pos, '\0', 0.0, search_radius);

            // no atoms are near here, continue to the next donor atom
            if (near_atoms.empty()) { std::cout << "No atoms nearby\n"; continue; }

            double min_distance = UINT16_MAX;
            gemmi::NeighborSearch::Mark *min_atom;

            for (const auto &atom: near_atoms) {

                // skip if the near atom is on the same residue
                if (atom->chain_idx == glycosite.chain_idx && atom->residue_idx == glycosite.residue_idx) {
                    continue;
                }
                gemmi::Residue bound_residue = m_structure->models[glycosite.model_idx].chains[atom->chain_idx].residues[atom->residue_idx];

                // skip if the near atom is on the same seqid (unlikely)
                // if (bound_residue.seqid == residue.seqid) {  continue; }

                // skip if the near atom is part of a unknown residue
                if (m_database.find(bound_residue.name) == m_database.end()) { continue; }

                if (bound_residue.atoms.empty()) { continue;}
                gemmi::Atom near_atom = bound_residue.atoms[atom->atom_idx];
                double distance = (donor_atom.pos - near_atom.pos).length();

                if (distance < min_distance) {
                    min_atom = atom;
                    min_distance = distance;
                }
            }

            // check we have changed the min_distance, if not, no suitable atoms were found
            if (min_distance == UINT16_MAX) { continue; }

            auto closest_site = Glycosite(*min_atom);

            gemmi::Residue closest_bound_residue = Sails::Utils::get_residue_from_glycosite(closest_site, m_structure);
            gemmi::Atom closest_atom = Sails::Utils::get_atom_from_glycosite(closest_site, m_structure);
            auto closest_residue_data = m_database[closest_bound_residue.name];
            auto acceptors = closest_residue_data.acceptors;

            // check if the closest atom is a known acceptor
            auto is_acceptor = [closest_atom](AtomSet& atom_set) { return atom_set.atom1 == closest_atom.name;};
            if (std::find_if(acceptors.begin(), acceptors.end(), is_acceptor) == acceptors.end()) {

                std::cout << "Closest atom to " << Utils::format_residue_key(&residue) << "-" << donor_atom.name <<
                    " is " << Utils::format_residue_key(&closest_bound_residue) << "-" << closest_atom.name << " which is not in ";
                for (const auto& a: acceptors) {
                    std::cout << a.atom1 << ",";
                }
                std::cout << std::endl;
                continue;
            };

            // Add the sugar, and then linkage
            // This is required to ensure the sugar objects live until the Glycan goes out of scope.
            glycan.add_sugar(closest_atom.name, closest_bound_residue.seqid.num.value, closest_site);

            // sugars are stored with keys which are the seqIds
            glycan.add_linkage(glycosite, closest_site);

            queue.push(closest_site);
        }
    }
}


Sails::Glycan Sails::Topology::find_glycan_topology(Glycosite &glycosite) {

    Glycan glycan = {m_structure, m_database, glycosite};

    std::queue<Glycosite> to_check({glycosite});

    while (!to_check.empty()) {
        auto current_site = to_check.front();
        to_check.pop();

        auto a = Utils::get_residue_from_glycosite(current_site, m_structure);
        find_residue_near_donor(current_site, glycan, to_check);
    }
    return glycan;
}
