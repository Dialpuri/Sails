//
// Created by Jordan Dialpuri on 12/07/2024.
//

#include "../include/sails-refine.h"


double Sails::TorsionAngleRefiner::calculate_penalty(double angle, double angle_mean, double angle_stddev,
                                                     double penalty_factor) {
    // int std_deviations_allowed = 2;
    // double range = std_deviations_allowed * angle_stddev;
    // double lower_bound = angle_mean - range;
    // double upper_bound = angle_mean + range;
    //
    // double deviation = 0;
    // if (angle < lower_bound) {
    //     deviation = lower_bound - angle;
    // } else {
    //     deviation = angle - upper_bound;
    // }
    //
    // double penalty = penalty_factor * pow(deviation, 2);
    // return penalty;
    //
    double angle_r = angle * M_PI / 180.0;
    double angle_mean_r = angle_mean * M_PI / 180.0;
    double angle_stddev_r = angle_stddev * M_PI / 180.0;
    double diff = angle_r - angle_mean_r;
    double delta = atan2(sin(diff), cos(diff)) ;
    double penalty = pow(delta, 2) / pow(angle_stddev_r,2);
    return penalty * penalty_factor;
}

double Sails::TorsionAngleRefiner::calculate_penalty_factor() const {
    switch (m_density->get_score_method()) {
        case atomwise:
            return 1e-2;
        case rscc:
            return 1e-5;
        default:
            return 0;
    }
}

double Sails::TorsionAngleRefiner::score_function(std::vector<double> &all_angles) {
    std::vector<double> angles = {all_angles[1], all_angles[2], all_angles[3]};
    std::vector<double> torsions = {all_angles[4], all_angles[5], all_angles[6]};

    gemmi::Residue residue = gemmi::Residue(m_reference_residue);
    gemmi::Transform superpose_result = Model::superpose_atoms(m_all_atoms, m_reference_atoms, all_angles[0], angles,
                                                               torsions);
    gemmi::transform_pos_and_adp(residue, superpose_result);
    SuperpositionResult result = {residue, superpose_result, m_reference_residue};

    const double score = -m_density->score_result(result);

    double penalty = 0;
    double penalty_factor = calculate_penalty_factor();
    for (int i = 0; i < 3; i++) {
        penalty += calculate_penalty(angles[i], m_angle_mean[i], m_angle_range[i], penalty_factor);
        penalty += calculate_penalty(torsions[i], m_torsion_mean[i], m_torsion_range[i], penalty_factor);
    }

    double bond_length_delta = std::abs(all_angles[0] - m_length);
    if (bond_length_delta > 0.3) {
        penalty += bond_length_delta * 1e5;
    }
    // std::cout << penalty << " " << score << " " << penalty_factor << std::endl;

    return score + penalty;
}

Sails::SuperpositionResult Sails::TorsionAngleRefiner::refine() {
    std::vector<double> initial_simplex = {
        m_length,
        m_angle_mean[0], m_angle_mean[1], m_angle_mean[2],
        m_torsion_mean[0], m_torsion_mean[1], m_torsion_mean[2]
    };

    // gemmi::Residue reference_residue = gemmi::Residue(m_reference_residue);
    // gemmi::Transform reference_superpose_result = Model::superpose_atoms(m_all_atoms, m_reference_atoms, m_length, m_angle_mean,
    //                                                            m_torsion_mean);
    // gemmi::transform_pos_and_adp(reference_residue, reference_superpose_result);
    // SuperpositionResult reference_result = {reference_residue, reference_superpose_result, m_reference_residue};
    //
    // const double initial_score = m_density->score_result(reference_result);
    // double penalty = 0;
    // for (int i = 0; i < 3; i++) {
    //     penalty += calculate_penalty(m_angle_mean[i], m_angle_mean[i], m_angle_range[i], calculate_penalty_factor());
    //     penalty += calculate_penalty(m_torsion_mean[i], m_torsion_mean[i], m_torsion_range[i], calculate_penalty_factor());
    // }

    auto lambda = [&](std::vector<double> &x) -> double {
        return this->score_function(x);
    };

    std::vector<double> final_simplex = nelder_mead::find_min(lambda, initial_simplex, true, {}, 1e-8, 1e-8, 100000,
                                                              100000);

    std::vector<double> final_angles = {
        final_simplex[1], final_simplex[2], final_simplex[3]
    };

    std::vector<double> final_torsions = {
        final_simplex[4], final_simplex[5], final_simplex[6]
    };

    gemmi::Residue residue = gemmi::Residue(m_reference_residue);
    gemmi::Transform final_result =
            Model::superpose_atoms(m_all_atoms, m_reference_atoms, final_simplex[0], final_angles, final_torsions);
    gemmi::transform_pos_and_adp(residue, final_result);
    SuperpositionResult result = {residue, final_result, m_reference_residue};

    // const double final_score = m_density->score_result(result);
    // double final_penalty = 0;
    // for (int i = 0; i < 3; i++) {
    //     final_penalty += calculate_penalty(final_angles[i], m_angle_mean[i], m_angle_range[i], calculate_penalty_factor());
    //     final_penalty += calculate_penalty(final_angles[i], m_torsion_mean[i], m_torsion_range[i], calculate_penalty_factor());
    // }
    //
    //
    // std::cout << std::endl <<  "Initial score: " << initial_score << " - penalty: " << penalty << std::endl;
    // std::vector<std::string> labels = {"length", "alpha", "beta", "gamma", "psi", "phi", "omega"};
    // std::cout << "\nLabel\tOriginal\tNew" << std::endl;
    // for (int i = 0; i < final_simplex.size(); i++) {
    //     std::cout << labels[i] << "\t" << initial_simplex[i] << "\t" << final_simplex[i] << std::endl;
    // }
    // std::cout << "Final score: " << final_score << " - penalty: " << final_penalty << std::endl;


    return result;
}
