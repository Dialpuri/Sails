//
// Created by Jordan Dialpuri on 12/07/2024.
//

#ifndef SAILS_REFINE_H
#define SAILS_REFINE_H

#include "sails-density.h"

#include "../third-party/nelder-mead.h"

#include <gemmi/modify.hpp>
#include <gemmi/model.hpp>

#include <chrono>

#include "sails-linkage.h"

namespace Sails {
    struct SuperpositionResult;

    struct TorsionAngleRefiner {
        TorsionAngleRefiner(
            std::vector<gemmi::Atom *> &atoms,
            std::vector<gemmi::Atom> &reference_atoms,
            Density& density,
            SuperpositionResult& superposition_result,
            double length,
            std::vector<double> & angle_mean,
            std::vector<double> & angle_range,
            std::vector<double> & torsion_mean,
            std::vector<double> & torsion_range

            ) : m_all_atoms(atoms),
            m_reference_atoms(reference_atoms),
            m_density(&density),
            m_reference_residue(superposition_result.reference_residue),
            m_length(length),
            m_angle_mean(angle_mean),
            m_torsion_mean(torsion_mean),
            m_angle_range(angle_range),
            m_torsion_range(torsion_range) {}


        double score_function(std::vector<double>& angles); // {alpha, beta, gamma, phi, psi, omega}

        SuperpositionResult refine();

    private:
        static double calculate_penalty(double angle, double angle_mean, double angle_stddev, double penalty_factor);

        double m_length;
        std::vector<double> m_angle_range;
        std::vector<double> m_torsion_range;
        std::vector<double> m_angle_mean;
        std::vector<double> m_torsion_mean;

        std::vector<gemmi::Atom *> m_all_atoms;
        std::vector<gemmi::Atom> m_reference_atoms;
        Density* m_density;
        gemmi::Residue m_reference_residue;
    };
};


#endif //SAILS_REFINE_H
