//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../../include/density/sails-density.h"
#include "../../include/sails-refine.h"
#include <clipper/contrib/edcalc.h>
#include <clipper/contrib/sfweight.h>
#include <clipper/minimol/minimol.h>

#include "src/include/sails-score.h"


double Sails::Density::score_residue(gemmi::Residue &residue, const DensityScoreMethod &method) {
    switch (method) {
        case atomwise:
            return atomwise_score(residue);
        case rscc:
            return rscc_score(residue);
        case rsr:
            return rsr_score(residue);
        // case dds:
        //     return check_difference_density(residue, TODO);
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
        case q:
            return q_score(result.new_residue);
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

// gemmi::Grid<> Sails::Density::calculate_density_for_box(gemmi::Residue &residue, gemmi::Box<gemmi::Position> &box) const {
//
//     gemmi::DensityCalculator<gemmi::IT92<float>, float> density_calculator;
//
//     gemmi::Position size = box.get_size();
//     gemmi::UnitCell dummy_cell = {size.x, size.y, size.z, 90, 90, 90};
//     density_calculator.grid.unit_cell = dummy_cell;
//     density_calculator.grid.nu = size.x;
//     density_calculator.grid.nv = size.y;
//     density_calculator.grid.nw = size.z;
//     density_calculator.grid.spacegroup = get_work_grid()->spacegroup;
//     density_calculator.grid.axis_order = get_work_grid()->axis_order;
//
//     density_calculator.d_min = 1;
//     density_calculator.initialize_grid();
//     for (auto &atom: residue.atoms) {
//         density_calculator.add_atom_density_to_grid(atom);
//     }
//     density_calculator.grid.symmetrize_sum();
//     return density_calculator.grid;
// }
//
// gemmi::Grid<> Sails::Density::calculate_density_for_grid(gemmi::Residue &residue) const {
//
//     gemmi::DensityCalculator<gemmi::C4322<float>, float> density_calculator;
//
//     density_calculator.grid.copy_metadata_from(*get_best_grid());
//     density_calculator.grid.spacing[0] = get_best_grid()->spacing[0];
//     density_calculator.grid.spacing[1] = get_best_grid()->spacing[1];
//     density_calculator.grid.spacing[2] = get_best_grid()->spacing[2];
//
//     density_calculator.d_min = get_resolution();
//     density_calculator.initialize_grid();
//     for (auto &atom: residue.atoms) {
//         density_calculator.add_atom_density_to_grid(atom);
//     }
//     density_calculator.grid.symmetrize_sum();
//     auto x =  density_calculator.grid;
//     return std::move(x);
// }
//
// gemmi::Grid<> Sails::Density::calculate_density_for_structure(gemmi::Structure &structure) const {
//     gemmi::DensityCalculator<gemmi::IT92<float>, float> density_calculator;
//
//     density_calculator.grid.copy_metadata_from(*get_best_grid());
//     density_calculator.grid.spacing[0] = get_best_grid()->spacing[0];
//     density_calculator.grid.spacing[1] = get_best_grid()->spacing[1];
//     density_calculator.grid.spacing[2] = get_best_grid()->spacing[2];
//
//     density_calculator.d_min = get_resolution();
//     density_calculator.initialize_grid();
//     density_calculator.add_model_density_to_grid(structure.models[0]);
//     density_calculator.grid.symmetrize_sum();
//     auto x =  density_calculator.grid;
//     return std::move(x);
// }

template <typename T>
T Sails::Density::calculate_rscc(std::vector<T> obs_values, std::vector<T> calc_values) {
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
template float Sails::Density::calculate_rscc<float>(std::vector<float> obs_values, std::vector<float> calc_values);
template double Sails::Density::calculate_rscc<double>(std::vector<double> obs_values, std::vector<double> calc_values);


float Sails::Density::rscc_score(gemmi::Residue &residue) const {
    if (residue.atoms.empty()) throw std::runtime_error("Residue is empty during RSCC check");

    gemmi::Box <gemmi::Position> box;
    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }
    box.add_margin(2);

    // gemmi::Grid<> calc = calculate_density_for_box(residue, box);
    gemmi::Grid<> calc = calculate_density_for_grid(residue);
    gemmi::Model model = Utils::create_model(residue);

    gemmi::NeighborSearch ns = {model, get_best_grid()->unit_cell, 2};
    ns.populate();
    // gemmi::Ccp4<> m;
    // m.grid = calc;
    // m.update_ccp4_header();
    // m.write_ccp4_map("calc.map");
    //
    // std::vector rs = {residue};
    // Utils::save_residues_to_file(rs, "res.pdb");

    const gemmi::Position min = box.minimum;
    const gemmi::Position max = box.maximum;

    std::vector<float> obs_values = {};
    std::vector<float> calc_values = {};

    for (double x = min.x; x <= max.x; x += get_best_grid()->spacing[0]) {
        for (double y = min.y; y <= max.y; y += get_best_grid()->spacing[1]) {
            for (double z = min.z; z <= max.z; z += get_best_grid()->spacing[2]) {
                gemmi::Position position = {x, y, z};
                auto nearest_atom = ns.find_atoms(position, '*', 0, 2);
                if (!nearest_atom.empty()) {
                    obs_values.emplace_back(get_best_grid()->interpolate_value(position));
                    calc_values.emplace_back(calc.interpolate_value(position));
                }
            }
        }
    }

    return calculate_rscc<float>(obs_values, calc_values);
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

    return calculate_rscc<float>(obs_values, calc_values);
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

int Sails::Density::check_difference_density(gemmi::Residue &residue, std::pair<float, float> map_stats) const {

    float threshold = map_stats.first - 2 * map_stats.second;

    std::set<std::string> ring_atoms = {
        "C1", "C2", "C3", "C4", "C5", "O5"
    };
    int i = 0;
    for (auto & atom : residue.atoms) {
        // if (ring_atoms.count(atom.name) == 0) continue;
        if (get_difference_grid()->interpolate_value(atom.pos) < threshold) {
            i++;
        }
    }
    return i;
    // gemmi::Box <gemmi::Position> box;
    // for (auto &atom: residue.atoms) {
    //     box.extend(atom.pos);
    // }
    //
    // const gemmi::Position max = box.maximum;
    // const gemmi::Position min = box.minimum;
    //
    // float sum = 0.0f;
    // int points = 0;
    // constexpr double step_size = 0.5;
    // for (double x = min.x; x <= max.x; x += step_size) {
    //     for (double y = min.y; y <= max.y; y += step_size) {
    //         for (double z = min.z; z <= max.z; z += step_size) {
    //             gemmi::Position position = {x, y, z};
    //             float value = get_difference_grid()->interpolate_value(position);
    //             sum += abs(value);
    //             points++;
    //         }
    //     }
    // }
    //
    // return sum / points;
}

float Sails::Density::score_atomic_position(const gemmi::Atom &atom) const {
    return score_position(atom.pos);
}


float Sails::Density::score_position(const gemmi::Position &pos) const {
    return get_work_grid()->interpolate_value(pos);
}

std::pair<float, float> Sails::Density::calculate_map_statistics(const gemmi::Grid<> *grid) const {
    const float sum = std::accumulate(grid->data.begin(), grid->data.end(), 0.0f);
    float mean = sum / grid->data.size();

    float sq_sum = std::accumulate(grid->data.begin(), grid->data.end(), 0.0,
            [mean](const double acc, const double x) {
                const double diff = x - mean;
                return acc + diff * diff;
            });

    float stdev = std::sqrt(sq_sum / grid->data.size());

    return std::make_pair(mean, stdev);
}

double Sails::Density::q_score(gemmi::Residue &residue) {
    auto [mean, stddev] = get_map_stats();

    const float A = mean + (10 * stddev);
    const float B = mean - stddev;
    constexpr float sigma = 0.6;
    constexpr int N = 8;

    gemmi::Model model = Utils::create_model(residue);
    gemmi::NeighborSearch ns = {model, get_best_grid()->unit_cell, 2};
    ns.populate();

    std::vector<double> residue_q_scores = {};

    for (int a = 0; a < residue.atoms.size(); a++) {
        Glycosite atom_site = {0, 0, 0, a};
        double atom_q = Score::QScore::calculate_q_score(residue.atoms[a].pos, atom_site, get_work_grid(),
            ns, A, B, sigma, N);
        residue_q_scores.emplace_back(atom_q);
    }

    const double mean_residue_q_score = std::accumulate(residue_q_scores.begin(), residue_q_scores.end(), 0.0)
                                        / static_cast<int>(residue.atoms.size());

    return mean_residue_q_score;
}
