//
// Created by Jordan Dialpuri on 08/08/2024.
//

#include "../include/sails-solvent.h"

gemmi::Box<gemmi::Position> Sails::SolventAccessibility::calculate_box() {
    gemmi::Box<gemmi::Position> structure_box;
    for (auto &m: m_structure->models) {
        for (auto &c: m.chains) {
            for (auto &r: c.residues) {
                for (auto &a: r.atoms) {
                    structure_box.extend(a.pos);
                }
            }
        }
    }
    structure_box.add_margin(8);
    return structure_box;
}

Sails::SolventAccessibility::SolventAccessibilityMap Sails::SolventAccessibility::create_solvent_accessibility_map() {
    SolventAccessibilityMap solvent_accessibility_map;
    for (int mi = 0; mi < m_structure->models.size(); mi++) {
        for (int ci = 0; ci < m_structure->models[mi].chains.size(); ci++) {
            for (int ri = 0; ri < m_structure->models[mi].chains[ci].residues.size(); ri++) {
                for (int ai = 0; ai < m_structure->models[mi].chains[ci].residues[ri].atoms.size(); ai++) {
                    Glycosite g = {mi, ci, ri, ai};
                    solvent_accessibility_map[g] = 0;
                }
            }
        }
    }
    return solvent_accessibility_map;
}

Sails::SolventAccessibility::SolventAccessibilityMap Sails::SolventAccessibility::average_solvent_accessibilty_map(
    const SolventAccessibilityMap &map) {
    SolventAccessibilityMap sum_map;
    SolventAccessibilityMap count_map;

    for (auto &[site, sa]: map) {
        Glycosite residue_site = {site.model_idx, site.chain_idx, site.residue_idx};
        sum_map[residue_site] += sa;
        count_map[residue_site] += 1;
    }

    SolventAccessibilityMap average_map;
    for (auto &[site, sum]: sum_map) {
        double count = count_map[site];
        if (count == 0 || sum == 0) {
            average_map[site] = 0;
            continue;
        }
        average_map[site] = sum / count;
    }
    return average_map;
}

std::map<Sails::Glycosite, double> Sails::SolventAccessibility::calculate_solvent_accessibility() {
    gemmi::Box<gemmi::Position> structure_box = calculate_box();
    auto min = structure_box.minimum;
    auto max = structure_box.maximum;
    double spacing = 1;

    double step_x = (max.x - min.x) / spacing;
    double step_y = (max.y - min.y) / spacing;
    double step_z = (max.z - min.z) / spacing;

    gemmi::NeighborSearch ns = {m_structure->models[0], m_structure->cell, 1.5};
    ns.populate();

    std::map<Glycosite, double> solvent_acessibility =create_solvent_accessibility_map();

    for (int x = 0; x < floor(step_x); x++) {
        for (int y = 0; y < floor(step_y); y++) {
            for (int z = 0; z < floor(step_z); z++) {
                gemmi::Position position = {min.x + x, min.y + y, min.z + z};
                double interpolated_value = m_solvent_mask.tricubic_interpolation(position);
                if (interpolated_value != 0 && interpolated_value != 8) {
                    auto near_atoms = ns.find_atoms(position, '\0', 0, 1.4);
                    for (const auto &near_atom: near_atoms) {
                        auto g = Glycosite(*near_atom);
                        // solvent_acessibility[g] += (8-interpolated_value);
                        solvent_acessibility[g] += 1;
                    }
                }
            }
        }
    }


    SolventAccessibilityMap flatmap = average_solvent_accessibilty_map(solvent_acessibility);
    return flatmap;
}

void Sails::SolventAccessibility::calculate_solvent_mask() {
    gemmi::SolventMasker masker = {gemmi::AtomicRadiiSet::Refmac, 1};
    constexpr int model_index = 0;

    masker.put_mask_on_grid(m_solvent_mask, m_structure->models[model_index]);
}
