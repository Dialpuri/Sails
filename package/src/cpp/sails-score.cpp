//
// Created by Jordan Dialpuri on 22/10/2025.
//


#include "../include/sails-score.h"

#include <gemmi/resinfo.hpp>

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

std::map<Sails::Glycosite, double> Sails::Score::calculate_qscores(Sails::Density *density, gemmi::Structure *structure,
    ResidueDatabase &residue_database) {

    constexpr double radius = 3;
    auto ns = gemmi::NeighborSearch(structure->models[0], structure->cell, radius);
    ns.populate();

    auto [mean, stddev] = density->calculate_map_statistics(density->get_best_grid());

    const float A = mean + (10 * stddev);
    const float B = mean - stddev;
    const float sigma = 0.6;
    constexpr int N = 8;

    std::map<Glycosite, double> qscores;
    for (auto & model : structure->models) {
        for (int c = 0; c < model.chains.size(); c++) {
            for (int r = 0; r < model.chains[c].residues.size(); r++) {
                const gemmi::Residue* residue_ptr = &model.chains[c].residues[r];
                if (residue_database.count(residue_ptr->name) > 0) {
                    const ResidueData& residue_data = residue_database.at(residue_ptr->name);
                    if (!residue_data.is_sugar) continue;
                    Glycosite site = {0, c, r, 0};

                    std::vector<double> residue_q_scores = {};
                    for (int a = 0; a < residue_ptr->atoms.size(); a++) {
                        Glycosite atom_site = {0, c, r, a};

                        double atom_q = Sails::Score::QScore::calculate_q_score(residue_ptr->atoms[a].pos, atom_site,
                            density->get_best_grid(), ns, A, B, sigma, N);
                        residue_q_scores.emplace_back(atom_q);
                    }
                    double mean_residue_q_score = std::accumulate(residue_q_scores.begin(), residue_q_scores.end(),
                        0.0) / static_cast<int>(residue_ptr->atoms.size());
                    qscores[site] = mean_residue_q_score;

                }
            }
        }
    }

    return qscores;
}

double Sails::Score::calculate_clash_score(gemmi::Residue *residue, gemmi::Structure *structure) {
    constexpr double radius = 1;
    gemmi::NeighborSearch ns = gemmi::NeighborSearch(structure->models[0], structure->cell, radius);
    ns.populate();

    double clash_score = 0;
    for (auto &atom: residue->atoms) {
        auto nearest_atoms = ns.find_atoms(atom.pos, '\0', 0, radius);
        clash_score += static_cast<double>(nearest_atoms.size());
    }
    return clash_score;
}

std::vector<gemmi::Position> Sails::Score::QScore::fibonacci_sphere(int samples, float radius, const gemmi::Position &center) {
    std::vector<gemmi::Position> positions;
    const double offset = 2.0 / samples;
    const double increment = M_PI * (3.0 - sqrt(5.0));

    for (int i = 0 ; i < samples; i++) {
        const double y = ((i * offset) - 1) + (offset / 2);
        const double r = sqrt(1 - pow(y,2));

        const double phi = i * increment;

        const double x = cos(phi) * r;
        const double z = sin(phi) * r;

        gemmi::Position position = {x, y, z};
        position *= radius;
        position += center;
        positions.emplace_back(position);
    }
    return positions;
}

std::vector<gemmi::Position> Sails::Score::QScore::get_radial_points(const gemmi::Position &position, float radius, int N,
                                                                     Glycosite &site, gemmi::NeighborSearch &ns) {

    std::vector<gemmi::Position> positions;
    constexpr int max_iter = 200;

    for (int i = 0 ; i < max_iter ; i++) {
        std::vector<gemmi::Position> sampled_sphere = fibonacci_sphere(N+i, radius, position);
        for (const auto& sampled_position: sampled_sphere) {
            // const gemmi::NeighborSearch::Mark* nearest_atom = ns.find_nearest_atom(sampled_position);
            // auto nearest_site = Glycosite(*nearest_atom);
            // if (nearest_site == site) {
            positions.emplace_back(sampled_position);
            // }

            if (positions.size() >= N) {
                break;
            }
        }
        if (positions.size() >= N) {
            break;
        }
    }
    return positions;
}

std::vector<double> Sails::Score::QScore::sample_density(const gemmi::Grid<> *grid, std::vector<gemmi::Position> &positions) {
    std::vector<double> values;
    for (auto& position: positions) {
        double value = grid->tricubic_interpolation(position);
        values.emplace_back(value);
    }
    return values;
}

double Sails::Score::QScore::calculate_q_score(const gemmi::Position & position, Glycosite &site, const gemmi::Grid<> *grid,
                                               gemmi::NeighborSearch &ns, float A, float B, float sigma, int N) {

    const int M = 21;
    std::vector<double> sample_space(M);
    for (int i = 0; i < M; i++)
        sample_space[i] = (2.0f / (M - 1)) * i;

    std::vector u(N, std::vector<double>(M, 0));
    std::vector v(N, std::vector<double>(M, 0));

    for (int i = 0; i < M; i++) {
        const double radius = sample_space[i];
        const double gaussian_sample = A * exp(-0.5 *  pow(radius / sigma, 2)) + B;

        auto radial_pts = get_radial_points(position, radius, N, site, ns);
        if (radial_pts.size() != static_cast<size_t>(N))
            continue;

        auto u_samples = sample_density(grid, radial_pts);

        for (int j = 0; j < N; j++) {
            u[j][i] = u_samples[j];
            v[j][i] = gaussian_sample;
        }
    }

    std::vector<double> u_flat, v_flat;
    for (int j = 0; j < N; j++) {
        const double mean_u = std::accumulate(u[j].begin(), u[j].end(), 0.0) / M;
        const double mean_v = std::accumulate(v[j].begin(), v[j].end(), 0.0) / M;
        for (int i = 0; i < M; i++) {
            u_flat.push_back(u[j][i] - mean_u);
            v_flat.push_back(v[j][i] - mean_v);
        }
    }

    double numerator = 0;
    double sum_u2 = 0;
    double sum_v2 = 0;

    for (size_t i = 0; i < u_flat.size(); i++) {
        numerator += u_flat[i] * v_flat[i];
        sum_u2 += u_flat[i] * u_flat[i];
        sum_v2 += v_flat[i] * v_flat[i];
    }

    return numerator / (std::sqrt(sum_u2) * std::sqrt(sum_v2));
}
