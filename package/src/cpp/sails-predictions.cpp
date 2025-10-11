//
// Created by Jordan Dialpuri on 07/10/2025.
//

#include "../include/sails-predictions.h"


std::optional<gemmi::NeighborSearch> Sails::Predictions::create_neighbour_search(
    gemmi::Grid<> *grid, float threshold, const gemmi::UnitCell &unit_cell) {

    gemmi::Model model = gemmi::Model(0);
    gemmi::Chain chain = gemmi::Chain("A");

    int seqid = 0;
    for (int u = 0; u < grid->nu; u++) {
        for (int v = 0; v < grid->nv; v++) {
            for (int w = 0; w < grid->nw; w++) {

                gemmi::Grid<>::Point point = grid->get_point(u, v, w);
                if (*point.value < threshold) {
                    continue;
                }
                gemmi::Position position = grid->point_to_position(point);
                gemmi::Atom atom;
                atom.name = "X";
                atom.element = gemmi::Element("C");
                atom.pos = position;
                gemmi::Residue residue = gemmi::Residue();
                residue.name = "PRD";
                residue.seqid = gemmi::SeqId(++seqid, '0');
                residue.atoms.emplace_back(atom);
                chain.residues.emplace_back(residue);
            }
        }
    }

    if (seqid == 0) {
        return std::nullopt;
    }

    model.chains = {chain};

    gemmi::NeighborSearch ns = {model, unit_cell, 2};
    ns.populate();
    return ns;
}

Sails::Glycosites Sails::Predictions::find_potential_sites(gemmi::Structure &structure) {
    if (m_glycan_map == nullptr) {
        throw std::invalid_argument("Glycan map is null");
    }
    if (m_protein_map == nullptr) {
        return find_potential_sites_using_glycan(structure);
    }
    return find_potential_sites_using_protein_glycan(structure);
}

Sails::Glycosites Sails::Predictions::find_potential_sites_using_glycan(gemmi::Structure &structure) {

    Glycosites potential_sites = {};

    std::optional<gemmi::NeighborSearch> ns_optional = create_neighbour_search(m_glycan_map, 0.1, structure.cell);
    if (!ns_optional.has_value()) {
        return potential_sites;
    }
    gemmi::NeighborSearch ns = ns_optional.value();

    for (int m = 0; m < structure.models.size(); m++) {
        for (int c = 0; c < structure.models[m].chains.size(); c++) {
            for (int r = 0; r < structure.models[m].chains[c].residues.size(); r++) {

                Glycosite site = {m, c, r};
                gemmi::Residue residue = structure.models[m].chains[c].residues[r];
                std::string residue_name = residue.name;

                if (protein_donors.find(residue_name) == protein_donors.end()) {
                    continue;
                }

                std::vector<AtomSet> donor_sets = m_residue_database[residue_name].donors;

                for (const auto& donor_set : donor_sets) {
                    std::string last_donor_atom_name = donor_set.atom3;
                    gemmi::Atom* last_donor_atom = residue.find_atom(last_donor_atom_name, '*');
                    auto nearby_points = ns.find_atoms(last_donor_atom->pos, '*', 0.1, 2);

                    if (nearby_points.empty()) {
                        continue;
                    }

                    potential_sites.emplace_back(site);
                    break;
                }
            }
        }
    }

    return potential_sites;
}

Sails::Glycosites Sails::Predictions::find_potential_sites_using_protein_glycan(gemmi::Structure &structure) {
    Glycosites potential_sites = {};

    std::optional<gemmi::NeighborSearch> ns_optional = create_neighbour_search(m_protein_map, 0.1, structure.cell);
    if (!ns_optional.has_value()) {
        return potential_sites;
    }
    gemmi::NeighborSearch ns = ns_optional.value();

    for (int m = 0; m < structure.models.size(); m++) {
        for (int c = 0; c < structure.models[m].chains.size(); c++) {
            for (int r = 0; r < structure.models[m].chains[c].residues.size(); r++) {

                Glycosite site = {m, c, r};
                gemmi::Residue residue = structure.models[m].chains[c].residues[r];
                std::string residue_name = residue.name;

                if (protein_donors.find(residue_name) == protein_donors.end()) {
                    continue;
                }

                std::vector<AtomSet> donor_sets = m_residue_database[residue_name].donors;

                for (const auto& donor_set : donor_sets) {
                    std::string last_donor_atom_name = donor_set.atom3;
                    gemmi::Atom* last_donor_atom = residue.find_atom(last_donor_atom_name, '*');
                    if (last_donor_atom == nullptr) {
                        continue;
                    }
                    auto nearby_points = ns.find_atoms(last_donor_atom->pos, '*', 0.1, 1);

                    if (nearby_points.empty()) {
                        continue;
                    }

                    potential_sites.emplace_back(site);
                    break;
                }
            }
        }
    }
    return potential_sites;
}
