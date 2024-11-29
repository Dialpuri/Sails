//
// Created by Jordan Dialpuri on 12/07/2024.
//

#ifndef SAILS_REFINE_H
#define SAILS_REFINE_H

#include "density/sails-density.h"
#include "sails-linkage.h"
#include "./third-party/nelder-mead.h"

#include <gemmi/modify.hpp>
#include <gemmi/model.hpp>

#include <chrono>


namespace Sails {
    struct SuperpositionResult;

    struct TorsionAngleRefiner {
        /**
         * Refines torsion angles for a set of atoms using reference atoms and density information.
         *
         * @param atoms A reference to a vector of Atom pointers.
         * @param reference_atoms A reference to a vector of Atom objects representing the reference atoms.
         * @param density A reference to a Density object.
         * @param superposition_result A reference to a SuperpositionResult object.
         * @param length A double representing the length.
         * @param angle_mean A reference to a vector of doubles representing the mean angles.
         * @param angle_range A reference to a vector of doubles representing the angle ranges.
         * @param torsion_mean A reference to a vector of doubles representing the mean torsions.
         * @param torsion_range A reference to a vector of doubles representing the torsion ranges.
         */
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


        /**
         * Calculates the score function for refining torsion angles based on a set of angles and the density score.
         *
         * @param all_angles A reference to a vector of doubles representing the input angles.
         *
         * Angles are in the order alpha beta gamma psi phi omega
         *
         * @return The calculated score function value.
         */
        double score_function(std::vector<double>& angles);

        /**
         * Refines torsion angles using a Nelder-Mead optimization algorithm.
         *
         * @return A SuperpositionResult object containing the refined residue, transform, and reference residue.
         */
        SuperpositionResult refine();



    private:

        /**
        * Calculates the penalty for a given angle based on its deviation from the mean angle
        * using the penalty factor.
        *
        * @param angle The angle to calculate the penalty for.
        * @param angle_mean The mean angle.
        * @param angle_stddev The standard deviation of the angle.
        * @param penalty_factor The penalty factor to be applied.
        *
        * @return The penalty value calculated for the given angle.
        */
        static double calculate_penalty(double angle, double angle_mean, double angle_stddev, double penalty_factor);

        /**
         * Calculates the penalty factor after looking at the scoring method in the m_density object.
         *
         * @return A double representing the penalty factor for the given density method.
         */
        [[nodiscard]] double calculate_penalty_factor() const;

        /**
         * The length of a linkage bond
         */
        double m_length;

        /**
         * Allowed angle range for alpha, beta and gamma linkage angles.
         */
        std::vector<double> m_angle_range;

        /**
         * Allowed angle range for phi, psi and omega linkage torsions.
         */
        std::vector<double> m_torsion_range;

        /**
         * Angle mean for alpha, beta and gamma linkage angles.
         */
        std::vector<double> m_angle_mean;

        /**
         * Torsion angle mean range for phi, psi and omega linkage torsions.
         */
        std::vector<double> m_torsion_mean;

        /**
         * List of pointers to all the atoms in the linkage to be refined
         */
        std::vector<gemmi::Atom *> m_all_atoms;

        /**
         * List of pointers to reference atoms in the linkage to be refined
         */
        std::vector<gemmi::Atom> m_reference_atoms;

        /**
         * Density object used for scoring to any derived density method (Xtal or EM)
         */
        Density* m_density;

        /**
         * The reference residue used for superpositons
         */
        gemmi::Residue m_reference_residue;
    };
};


#endif //SAILS_REFINE_H
