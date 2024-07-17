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
            double length
            ) : m_all_atoms(atoms),
            m_reference_atoms(reference_atoms),
            m_density(&density),
            m_reference_residue(superposition_result.reference_residue),
            m_length(length){}


        double score_function(std::vector<double>& angles); // {alpha, beta, gamma, phi, psi, omega}

        SuperpositionResult refine(std::vector<double> &initial_angles, std::vector<double> &initial_torsions);

    private:
        double m_length;
        std::vector<gemmi::Atom *> m_all_atoms;
        std::vector<gemmi::Atom> m_reference_atoms;
        Density* m_density;
        gemmi::Residue m_reference_residue;
    };
};


#endif //SAILS_REFINE_H
