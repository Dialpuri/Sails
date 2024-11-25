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

    double mean_delta = abs(angle - angle_mean);
    if (mean_delta < range) {
        return penalty_factor * mean_delta;
    }

    double deviation = 0;
    if (angle < lower_bound) {
        deviation = lower_bound - angle;
    } else {
        deviation = angle - upper_bound;
    }
    return deviation;
}

double Sails::TorsionAngleRefiner::calculate_penalty_factor() const {
    switch (m_density->get_score_method()) {
        case atomwise:
            return 1e-3;
        case rscc:
            return 1e-5;
        default:
            return 0;
    }
}

double Sails::TorsionAngleRefiner::score_function(std::vector<double> &all_angles) {
    std::vector<double> angles = {all_angles[0], all_angles[1], all_angles[2]};
    std::vector<double> torsions = {all_angles[3], all_angles[4], all_angles[5]};

    gemmi::Residue residue = gemmi::Residue(m_reference_residue);
    gemmi::Transform superpose_result = Model::superpose_atoms(m_all_atoms, m_reference_atoms, m_length, angles,
                                                               torsions);
    gemmi::transform_pos_and_adp(residue, superpose_result);
    SuperpositionResult result = {residue, superpose_result, m_reference_residue};

    const double score = m_density->score_result(result);

    double penalty = 0;
    double penalty_factor = calculate_penalty_factor();
    for (int i = 0; i < 3; i++) {
        penalty += calculate_penalty(angles[i], m_angle_mean[i], m_angle_range[i], penalty_factor);
        penalty += calculate_penalty(torsions[i], m_torsion_mean[i], m_torsion_range[i], penalty_factor);
    }

    return penalty - score;
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


std::vector<double> Sails::GlycanRegulariser::simplex(std::function<double(std::vector<double> &)> function,
                                                      std::vector<std::vector<double> > &arguments, int no_parameters,
                                                      int max_cycles) {
    enum STEP { UNKN, EXTN, NRML, CTRN, CTRX };
    double tolerance = 1e-8;

    std::vector<std::vector<double> > parameters(no_parameters + 1);
    std::vector<double> results(no_parameters + 1);

    for (int i = 0; i < arguments.size(); ++i) {
        parameters[i] = arguments[i];
        results[i] = function(arguments[i]);
    }

    int worst_index = 0;
    int best_index = 0;

    double f0 = 0.0;
    double f1 = 0.0;
    double f2 = 0.0;

    for (int cycle_index = 0; cycle_index < max_cycles; ++cycle_index) {
        worst_index = 0;
        best_index = 0;

        for (int i = 0; i < parameters.size(); ++i) {
            if (results[i] > results[worst_index]) { worst_index = i; }
            if (results[i] < results[best_index]) { best_index = i; }
        }

        if (results[worst_index] - results[best_index] < tolerance) {
            // std::cout << "Results[worst_index]: " << results[worst_index] << "Results[best_index]: " << results[
            //             best_index]
            //         << " diff " << results[worst_index] - results[best_index] << std::endl;
            break;
        }

        std::vector<double> centroid(no_parameters, 0.0);
        std::vector<double> shift(no_parameters);
        std::vector<double> t0(no_parameters);
        std::vector<double> t1(no_parameters);
        std::vector<double> t2(no_parameters);
        double weight_sum = 0.0;

        for (int i = 0; i < no_parameters; ++i) {
            if (i == worst_index) continue;

            double weight = 1.0;
            for (int j = 0; j < no_parameters; ++j) {
                centroid[j] = weight * parameters[i][j];
            }
            weight_sum += weight;
        }

        for (int j = 0; j < no_parameters; ++j) {
            centroid[j] /= weight_sum;
            shift[j] = centroid[j] - parameters[worst_index][j];
        }

        for (int j = 0; j < no_parameters; ++j) {
            t0[j] = centroid[j] - 0.5 * shift[j];
            t1[j] = centroid[j] + 1.0 * shift[j];
            t2[j] = centroid[j] + 2.0 * shift[j];
        }

        f1 = function(t1);
        STEP step = UNKN;

        if (f1 < results[worst_index]) {
            if (f1 < results[best_index]) {
                f2 = function(t2);
                if (f2 < f1) {
                    step = EXTN;
                } else {
                    step = NRML;
                }
            } else {
                step = NRML;
            }
        }

        if (step == UNKN) {
            f0 = function(t0);
            if (f0 < results[worst_index]) {
                step = CTRN;
            } else {
                step = CTRX;
            }
        }

        if (step == EXTN) {
            parameters[worst_index] = t2;
            results[worst_index] = f2;
        } else if (step == NRML) {
            parameters[worst_index] = t1;
            results[worst_index] = f1;
        } else if (step == CTRN) {
            parameters[worst_index] = t0;
            results[worst_index] = f0;
        } else {
            for (int i = 0; i < parameters.size(); ++i) {
                if (i == best_index) continue;
                for (int j = 0; j < no_parameters; ++j) {
                    parameters[i][j] = 0.5 * (parameters[i][j] + parameters[best_index][j]);
                }
                results[i] = function(parameters[i]);
            }
        }
    }
    return parameters[best_index];
}


void Sails::GlycanRegulariser::regularise(Glycan &glycan) {
    std::vector<LinkageData> linkage_data = glycan.extract_linkage_data(m_linkage_database);

    Linkage linkage = glycan.linkage_list[0];

    int parameters = glycan.linkage_list.size() * 2;
    std::vector<double> initial_values = {};

    for (int i = 0; i < glycan.linkage_list.size(); i++) {
        auto priority_cluster = std::find_if(linkage_data[i].clusters.begin(), linkage_data[i].clusters.end(),
                                             [](const Cluster &cluster) {
                                                 return cluster.priority;
                                             });
        initial_values.emplace_back(priority_cluster->torsions.psi.mean);
        initial_values.emplace_back(priority_cluster->torsions.phi.mean);
    }

    std::vector<std::vector<double> > initial_simplex = {initial_values};
    initial_simplex.reserve(parameters);
    for (int i = 0; i < parameters; i++) {
        int step = 2;
        std::vector<double> current_simplex(initial_values.begin(), initial_values.end());
        current_simplex[i] += step;
        initial_simplex.emplace_back(current_simplex);
    }

    gemmi::NeighborSearch ns = glycan.create_amino_acid_neighbor_search();

    auto lambda = [&](std::vector<double> &x) -> double {
        gemmi::Structure trial_structure(glycan.get_structure());
        double penalty = 0.0;
        int value_index = -1;

        for (int i = 0; i < glycan.linkage_list.size(); i++) {

            const auto priority_cluster = std::find_if(linkage_data[i].clusters.begin(), linkage_data[i].clusters.end(),
                                                 [](const Cluster &cluster) {return cluster.priority;});

            TorsionSet trial_torsion = priority_cluster->torsions;
            trial_torsion.psi.mean = x[++value_index];
            trial_torsion.phi.mean = x[++value_index];
            Cluster angles = {
                priority_cluster->angles,
                trial_torsion,
                true
            };

            Model::update_linkage_torsion(glycan, glycan.linkage_list[i], angles, &trial_structure,
                                          m_residue_database);

            std::vector<double> database_means = priority_cluster->torsions.get_means_in_order();
            std::vector<double> database_stddevs = priority_cluster->torsions.get_stddev_in_order();
            std::vector<double> trial_means = trial_torsion.get_means_in_order();;
            for (int j = 0; j < 2; ++j) {
                penalty += calculate_penalty(trial_means[j], database_means[j], database_stddevs[j], 1e-2);
            }

        }

        double clash_score = Model::calculate_glycan_clash_score(glycan, &ns, &trial_structure);
        return clash_score + penalty;
    };

    std::vector<double> final_values = simplex(lambda, initial_simplex, parameters, 10000);
    // std::vector<double> final_values = nelder_mead::find_min(lambda, initial_values, true, {}, 1e-8, 1e-8,
    // 100000,
    // 100000);

    int value_index = -1;
    gemmi::Structure trial_structure(glycan.get_structure());

    for (int i = 0; i < glycan.linkage_list.size(); i++) {
        TorsionSet t = linkage_data[i].clusters[0].torsions;
        t.psi.mean = final_values[++value_index];
        t.phi.mean = final_values[++value_index];
        Cluster angles = {
            linkage_data[i].clusters[0].angles,
            t,
            false
        };
        Model::update_linkage_torsion(glycan, glycan.linkage_list[i], angles, &trial_structure,
                                      m_residue_database);
    }


    // exit(-1);
}
