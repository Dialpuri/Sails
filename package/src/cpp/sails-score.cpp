//
// Created by Jordan Dialpuri on 22/10/2025.
//


#include "../include/sails-score.h"

#include "src/include/sails-utils.h"

std::map<Sails::Glycosite, double> Sails::Score::calculate_rsccs(Density *density, gemmi::Structure *structure, ResidueDatabase &residue_database) {
    gemmi::Grid<> calculated_density = density->calculate_density_for_structure(*structure);

    constexpr double radius = 2;
    auto ns = gemmi::NeighborSearch(structure->models[0], structure->cell, radius);
    ns.populate();

    gemmi::Grid<> best_grid = *density->get_best_grid();
    std::map<Glycosite, std::vector<std::pair<double, double>>> residue_pairs;

    for (auto point: best_grid) {
        gemmi::Position  position = best_grid.point_to_position(point);
        auto mark = ns.find_nearest_atom(position, radius);
        if (mark == nullptr) continue;

        auto site = Glycosite(0, mark->chain_idx, mark->residue_idx, 0);
        const gemmi::Residue* residue_ptr = &structure->models[site.model_idx].chains[site.chain_idx].residues[site.residue_idx];
        if (residue_database.count(residue_ptr->name) == 0) continue;
        const ResidueData& residue = residue_database.at(residue_ptr->name);
        if (!residue.is_sugar) continue;

        double obs = *point.value;
        double calc = calculated_density.interpolate_value(position);
        residue_pairs[site].emplace_back(obs, calc);

    }
    std::map<Glycosite, double> rsccs;

    for (const auto& [site, data]: residue_pairs) {
        auto [obs_values, calc_values] = Sails::Utils::split_pairs<double>(data);
        if (obs_values.empty() || calc_values.empty()) continue;

        rsccs[site] = Sails::Density::calculate_rscc<double>(obs_values, calc_values);
    }
    return rsccs;
}
