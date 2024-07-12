//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../include/sails-density.h"

#include <gemmi/calculate.hpp>
#include <src/include/sails-utils.h>
#include <gemmi/ccp4.hpp>


Sails::Density::Density(const std::string &mtz_path, const std::string& f_col, const std::string& phi_col) {
    gemmi::Mtz mtz = gemmi::read_mtz_file(mtz_path);
    m_mtz = std::move(mtz);
    m_grid = load_grid(m_mtz, f_col, phi_col, false);
    initialise_density_calculator();
}

Sails::Density::Density(gemmi::Mtz &mtz, const std::string& f_col, const std::string& phi_col) {
    m_mtz = std::move(mtz);
    m_grid = load_grid(m_mtz, f_col, phi_col, false);
    initialise_density_calculator();
}

double Sails::Density::score_residue(gemmi::Residue &residue, const DensityScoreMethod &method) {
    switch (method) {
        case atomwise:
            return atomwise_score(residue);
        case rscc:
            return rscc_score(residue);
        case rsr:
            return rsr_score(residue);
        default:
            return -1;
    }
}

gemmi::Grid<> Sails::Density::load_grid(const gemmi::Mtz &mtz, const std::string &f_col, const std::string &phi_col,
                                        bool normalise) {
    constexpr std::array<int, 3> null_size = {0, 0, 0};
    constexpr double sample_rate = 0;
    constexpr auto order = gemmi::AxisOrder::XYZ;

    const gemmi::Mtz::Column &f = mtz.get_column_with_label(f_col);
    const gemmi::Mtz::Column &phi = mtz.get_column_with_label(phi_col);
    const gemmi::FPhiProxy fphi(gemmi::MtzDataProxy{mtz}, f.idx, phi.idx);
    gemmi::Grid<> grid = gemmi::transform_f_phi_to_map2<float>(fphi, null_size, sample_rate, null_size, order);
    if (normalise) grid.normalize();

    return grid;
}

float Sails::Density::atomwise_score(const gemmi::Residue &residue) const {
    return std::transform_reduce(residue.atoms.begin(), residue.atoms.end(), 0.0f, std::plus<>(),
                                 [&](const gemmi::Atom &current_atom) {
                                     return m_grid.interpolate_value(current_atom.pos);
                                 });
}

void Sails::Density::calculate_density_for_box(gemmi::Residue &residue, gemmi::Box<gemmi::Position> &box) {
    density_calculator.initialize_grid();
    for (auto& atom: residue.atoms) {
        density_calculator.add_atom_density_to_grid(atom);
        expand_box(atom, box);
    }
    density_calculator.grid.symmetrize_sum();

}

float Sails::Density::rscc_score(gemmi::Residue &residue) {
    gemmi::Box<gemmi::Position> box;
    calculate_density_for_box(residue, box);

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    std::vector<float> obs_values = {};
    std::vector<float> calc_values = {};

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                obs_values.emplace_back(m_grid.interpolate_value(position));
                calc_values.emplace_back(density_calculator.grid.interpolate_value(position));
            }
        }
    }

    if (obs_values.size() != calc_values.size()) throw std::runtime_error("RSCC obs and calc lists are different sizes");

    float obs_average = std::accumulate(obs_values.begin(), obs_values.end(), 0.0f) / obs_values.size();
    float calc_average = std::accumulate(calc_values.begin(), calc_values.end(), 0.0f) / calc_values.size();

    float numerator = 0.0f;
    float obs_sum_sq = 0.0f;
    float calc_sum_sq = 0.0f;

    for (int i = 0; i < obs_values.size(); i++) {
        float obs_delta = obs_values[i]-obs_average;
        float calc_delta = calc_values[i]-calc_average;

        numerator += obs_delta*calc_delta;
        obs_sum_sq += (obs_delta*obs_delta);
        calc_sum_sq += (calc_delta*calc_delta);
    }

    float denominator = sqrt(obs_sum_sq*calc_sum_sq);

    if (denominator == 0.0f) throw std::runtime_error("RSCC Denominator is 0");
    return numerator / denominator;
}

float Sails::Density::rsr_score(gemmi::Residue &residue) {
    gemmi::Box<gemmi::Position> box;
    calculate_density_for_box(residue, box);

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    float numerator = 0.0f;
    float denominator = 0.0f;

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                float obs = m_grid.interpolate_value(position); // obs
                float calc = density_calculator.grid.interpolate_value(position); //calc
                numerator += abs(obs - calc);
                denominator += abs(obs + calc);
            }
        }
    }

    if (denominator == 0.0f) throw std::runtime_error("Box is empty");
    return numerator / denominator;
}


