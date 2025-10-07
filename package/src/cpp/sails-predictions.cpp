//
// Created by Jordan Dialpuri on 07/10/2025.
//

#include "../include/sails-predictions.h"


gemmi::NeighborSearch Sails::Predictions::create_neighbour_search(float threshold, gemmi::UnitCell unit_cell) {

    gemmi::Model model = gemmi::Model(0);
    gemmi::Chain chain = gemmi::Chain("A");

    int seqid = 0;
    for (int u = 0; u < m_glycan_map.nu; u++) {
        for (int v = 0; v < m_glycan_map.nv; v++) {
            for (int w = 0; w < m_glycan_map.nw; w++) {

                gemmi::Grid<>::Point point = m_glycan_map.get_point(u, v, w);
                if (*point.value < threshold) {
                    continue;
                }
                gemmi::Position position = m_glycan_map.point_to_position(point);
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

    model.chains = {chain};

    gemmi::Structure s;
    s.cell = m_glycan_map.unit_cell;
    s.spacegroup_hm = m_glycan_map.spacegroup->hm;
    s.models = {model};
    Utils::save_structure_to_file(s, "points.cif");

    std::cout << m_glycan_map.unit_cell.a << " " << m_glycan_map.unit_cell.b << " " << m_glycan_map.unit_cell.c << " "
    << m_glycan_map.unit_cell.alpha << " " << m_glycan_map.unit_cell.beta << " " << m_glycan_map.unit_cell.gamma << std::endl;
    gemmi::NeighborSearch ns = {model, unit_cell, 2};
    ns.populate();
    return ns;
}

Sails::Glycosites Sails::Predictions::find_potential_sites(gemmi::Structure &structure) {

    Glycosites potential_sites = {};

    gemmi::NeighborSearch ns = create_neighbour_search(0.1, structure.cell);

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
