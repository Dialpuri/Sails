//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../../include/density/sails-density.h"
#include "../../include/sails-refine.h"
#include <clipper/contrib/edcalc.h>
#include <clipper/contrib/sfweight.h>
#include <clipper/minimol/minimol.h>


double Sails::Density::score_residue(gemmi::Residue &residue, const DensityScoreMethod &method) {
    switch (method) {
        case atomwise:
            return atomwise_score(residue);
        case rscc:
            return rscc_score(residue);
        case rsr:
            return rsr_score(residue);
        case dds:
            return difference_density_score(residue);
        default:
            return -1;
    }
}

double Sails::Density::score_result(SuperpositionResult& result) {
    switch (get_score_method()) {
        case atomwise:
            return atomwise_score(result.new_residue);
        case rscc:
            return rscc_score(result);
        case rsr:
            return rsr_score(result);
        case dds:
            return difference_density_score(result.new_residue);
        default:
            return -1;
    }
}

float Sails::Density::atomwise_score(const gemmi::Residue &residue) const {
    return std::transform_reduce(residue.atoms.begin(), residue.atoms.end(), 0.0f, std::plus<>(),
                                 [&](const gemmi::Atom &current_atom) {
                                     return get_work_grid()->interpolate_value(current_atom.pos);
                                 }) / (residue.atoms.size());
}

gemmi::Grid<> Sails::Density::calculate_density_for_box(gemmi::Residue &residue, gemmi::Box<gemmi::Position> &box) const {

    gemmi::DensityCalculator<gemmi::IT92<float>, float> density_calculator;

    gemmi::Position size = box.get_size();
    gemmi::UnitCell dummy_cell = {size.x, size.y, size.z, 90, 90, 90};
    density_calculator.grid.unit_cell = dummy_cell;
    density_calculator.grid.nu = size.x;
    density_calculator.grid.nv = size.y;
    density_calculator.grid.nw = size.z;
    density_calculator.grid.spacegroup = get_work_grid()->spacegroup;
    density_calculator.grid.axis_order = get_work_grid()->axis_order;

    density_calculator.d_min = 1;
    density_calculator.initialize_grid();
    for (auto &atom: residue.atoms) {
        density_calculator.add_atom_density_to_grid(atom);
    }
    density_calculator.grid.symmetrize_sum();
    return density_calculator.grid;
}

gemmi::Grid<> Sails::Density::calculate_density_for_grid(gemmi::Residue &residue) const {

    gemmi::DensityCalculator<gemmi::C4322<float>, float> density_calculator;

    density_calculator.grid.copy_metadata_from(*get_work_grid());
    density_calculator.grid.spacing[0] = get_work_grid()->spacing[0];
    density_calculator.grid.spacing[1] = get_work_grid()->spacing[1];
    density_calculator.grid.spacing[2] = get_work_grid()->spacing[2];

    density_calculator.d_min = get_resolution();
    density_calculator.initialize_grid();
    for (auto &atom: residue.atoms) {
        density_calculator.add_atom_density_to_grid(atom);
    }
    density_calculator.grid.symmetrize_sum();
    auto x =  density_calculator.grid;
    return std::move(x);
}

float Sails::Density::calculate_rscc(std::vector<float> obs_values, std::vector<float> calc_values) {
    if (obs_values.size() != calc_values.size())
        throw std::runtime_error("RSCC obs and calc lists are different sizes");

    if (obs_values.empty()) throw std::runtime_error("Observation list is empty");
    if (calc_values.empty()) throw std::runtime_error("Calculated list is empty");

    float obs_average = std::accumulate(obs_values.begin(), obs_values.end(), 0.0f) / obs_values.size();
    float calc_average = std::accumulate(calc_values.begin(), calc_values.end(), 0.0f) / calc_values.size();

    if (calc_average == 0.0f) throw std::runtime_error("Calculated map average is 0");

    float numerator = 0.0f;
    float obs_sum_sq = 0.0f;
    float calc_sum_sq = 0.0f;

    for (int i = 0; i < obs_values.size(); i++) {
        float obs_delta = obs_values[i] - obs_average;
        float calc_delta = calc_values[i] - calc_average;

        numerator += obs_delta * calc_delta;
        obs_sum_sq += (obs_delta * obs_delta);
        calc_sum_sq += (calc_delta * calc_delta);
    }

    float denominator = sqrt(obs_sum_sq * calc_sum_sq);

    if (denominator == 0.0f) throw std::runtime_error("RSCC Denominator is 0");
    return numerator / denominator;
}


float Sails::Density::rscc_score(gemmi::Residue &residue) const {
    if (residue.atoms.empty()) throw std::runtime_error("Residue is empty during RSCC check");

    gemmi::Box <gemmi::Position> box;
    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }
    box.add_margin(1);

    // gemmi::Grid<> calc = calculate_density_for_box(residue, box);
    gemmi::Grid<> calc = calculate_density_for_grid(residue);

    // gemmi::Ccp4<> m;
    // m.grid = calc;
    // m.update_ccp4_header();
    // m.write_ccp4_map("calc.map");
    //
    // std::vector rs = {residue};
    // Utils::save_residues_to_file(rs, "res.pdb");


    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    std::vector<float> obs_values = {};
    std::vector<float> calc_values = {};

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                obs_values.emplace_back(get_best_grid()->interpolate_value(position));
                calc_values.emplace_back(calc.interpolate_value(position));
            }
        }
    }

    return calculate_rscc(obs_values, calc_values);
}

float Sails::Density::rscc_score(SuperpositionResult &result) {
    gemmi::Box <gemmi::Position> box;
    gemmi::Residue residue = result.new_residue;

    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }
    box.add_margin(1);

    auto calculated_maps = get_calculated_maps();
    if (calculated_maps->find(residue.name) == calculated_maps->end()) {
        gemmi::Grid<> reference = calculate_density_for_grid(result.reference_residue);
        calculated_maps->operator[](residue.name) = std::move(reference);
    }
    gemmi::Grid<> *calculated = &calculated_maps->operator[](residue.name);

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    std::vector<float> obs_values = {};
    std::vector<float> calc_values = {};

    gemmi::Residue r1, r2, r3;

    constexpr double step_size = 1;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                obs_values.emplace_back(get_best_grid()->interpolate_value(position));
                gemmi::Vec3 translated_position = result.transformation.inverse().apply(position);
                calc_values.emplace_back(calculated->interpolate_value(gemmi::Position(translated_position)));
            }
        }
    }

    return calculate_rscc(obs_values, calc_values);
}

float Sails::Density::rsr_score(gemmi::Residue &residue) {
    gemmi::Box <gemmi::Position> box;
    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }

    gemmi::Grid<> calc_grid = calculate_density_for_box(residue, box);

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    float numerator = 0.0f;
    float denominator = 0.0f;

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                float obs = get_work_grid()->interpolate_value(position);
                float calc = calc_grid.interpolate_value(position);
                numerator += abs(obs - calc);
                denominator += abs(obs + calc);
            }
        }
    }

    if (denominator == 0.0f) throw std::runtime_error("Box is empty");
    return numerator / denominator;
}

float Sails::Density::rsr_score(SuperpositionResult &result) {
    gemmi::Box <gemmi::Position> box;
    gemmi::Residue residue = result.new_residue;

    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }
    box.add_margin(1);

    // calculate map if not found
    auto calculated_maps = get_calculated_maps();
    if (calculated_maps->find(residue.name) == calculated_maps->end()) {
        gemmi::Grid<> reference = calculate_density_for_grid(result.reference_residue);
        calculated_maps->operator[](residue.name) = std::move(reference);
    }
    gemmi::Grid<> *calculated = &calculated_maps->operator[](residue.name);

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    float numerator = 0.0f;
    float denominator = 0.0f;

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                float obs = get_work_grid()->interpolate_value(position);
                gemmi::Vec3 translated_position = result.transformation.inverse().apply(position);
                float calc = calculated->interpolate_value(gemmi::Position(translated_position));
                numerator += abs(obs - calc);
                denominator += abs(obs + calc);
            }
        }
    }
    if (denominator == 0.0f) throw std::runtime_error("Box is empty");
    return numerator / denominator;
}

float Sails::Density::difference_density_score(gemmi::Residue &residue) const {
    gemmi::Box <gemmi::Position> box;
    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    float sum = 0.0f;
    int points = 0;
    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                float value = get_difference_grid()->interpolate_value(position);
                sum += abs(value);
                points++;
            }
        }
    }

    return sum / points;
}

float Sails::Density::score_atomic_position(const gemmi::Atom &atom) const {
    return score_position(atom.pos);
}


float Sails::Density::score_position(const gemmi::Position &pos) const {
    return get_work_grid()->interpolate_value(pos);
}
