//
// Created by Jordan Dialpuri on 12/07/2024.
//

#include "../include/sails-refine.h"

double Sails::TorsionAngleRefiner::score_function(std::vector<double> &all_angles) {
    std::vector<double> angles = {all_angles[0], all_angles[1], all_angles[2]};
    std::vector<double> torsions = {all_angles[3], all_angles[4], all_angles[5]};

    gemmi::Residue residue = gemmi::Residue(m_reference_residue);

    gemmi::Transform superpose_result = Model::superpose_atoms(m_all_atoms, m_reference_atoms, m_length, angles,
                                                               torsions);
    gemmi::transform_pos_and_adp(residue, superpose_result);
    SuperpositionResult result = {residue, superpose_result, m_reference_residue};
    // return -m_density.rscc_score(result);
    return -m_density->atomwise_score(residue);
}

Sails::SuperpositionResult Sails::TorsionAngleRefiner::refine(
    std::vector<double> &initial_angles,
    std::vector<double> &initial_torsions) {

    std::vector<double> initial_simplex = {
        initial_angles[0], initial_angles[1], initial_angles[2],
        initial_torsions[0], initial_torsions[1], initial_torsions[2]
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
    return result;
}
