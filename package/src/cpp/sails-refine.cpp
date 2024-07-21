//
// Created by Jordan Dialpuri on 12/07/2024.
//

#include "../include/sails-refine.h"


double Sails::TorsionAngleRefiner::calculate_penalty(double angle, double angle_mean, double angle_stddev,
                                                     double penalty_factor) {
    int std_deviations_allowed = 1;
    double range = std_deviations_allowed * angle_stddev;
    double lower_bound = angle_mean - range;
    double upper_bound = angle_mean + range;

    double deviation = 0;
    if (angle < lower_bound) {
        deviation = lower_bound - angle;
    } else {
        deviation = angle - upper_bound;
    }

    double penalty = penalty_factor * pow(deviation, 2);
    return penalty;
}

double Sails::TorsionAngleRefiner::score_function(std::vector<double> &all_angles) {
    std::vector<double> angles = {all_angles[0], all_angles[1], all_angles[2]};
    std::vector<double> torsions = {all_angles[3], all_angles[4], all_angles[5]};

    gemmi::Residue residue = gemmi::Residue(m_reference_residue);
    gemmi::Transform superpose_result = Model::superpose_atoms(m_all_atoms, m_reference_atoms, m_length, angles,
                                                               torsions);
    gemmi::transform_pos_and_adp(residue, superpose_result);
    SuperpositionResult result = {residue, superpose_result, m_reference_residue};

    double score = m_density->atomwise_score(residue);

    double penalty = 0;
    for (int i = 0; i < 3; i++) {
        penalty += calculate_penalty(angles[i], m_angle_mean[i], m_angle_range[i], 1e-3);
        penalty += calculate_penalty(torsions[i], m_torsion_mean[i], m_torsion_range[i], 1e-3);
    }

    return penalty-score;
}

Sails::SuperpositionResult Sails::TorsionAngleRefiner::refine() {
    std::vector<double> initial_simplex = {
        m_angle_mean[0], m_angle_mean[1], m_angle_mean[2],
        m_torsion_mean[0], m_torsion_mean[1], m_torsion_mean[2]
    };

    auto lambda = [&](std::vector<double> &x) -> double {
        return this->score_function(x);
    };

    std::vector<double> final_simplex = nelder_mead::find_min(lambda, initial_simplex, true, {}, 1e-8, 1e-8, 100000,
                                                              100000);

    std::vector<double> final_angles = {
        final_simplex[0], final_simplex[1], final_simplex[2]
    };

    std::vector<double> final_torsions = {
        final_simplex[3], final_simplex[4], final_simplex[5]
    };

    gemmi::Residue residue = gemmi::Residue(m_reference_residue);
    gemmi::Transform final_result =
            Model::superpose_atoms(m_all_atoms, m_reference_atoms, m_length, final_angles, final_torsions);
    gemmi::transform_pos_and_adp(residue, final_result);
    SuperpositionResult result = {residue, final_result, m_reference_residue};

    // std::vector<std::string> labels = {"alpha", "beta", "gamma", "psi", "phi", "omega"};
    // std::cout << "\nLabel\tOriginal\tNew" << std::endl;
    // for (int i = 0; i < final_simplex.size(); i++) {
    //     std::cout << labels[i] << "\t" << initial_simplex[i] << "\t" << final_simplex[i] << std::endl;
    // }

    return result;
}

