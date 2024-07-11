//
// Created by Jordan Dialpuri on 05/07/2024.
//

#ifndef SAILS_SAILS_VECTOR_H
#define SAILS_SAILS_VECTOR_H

#include <gemmi/unitcell.hpp>
#include "gemmi/math.hpp"

#include <clipper/core/coords.h>
#include <clipper/core/rotation.h>

namespace Sails {
    /**
     * Calculates the projected point given three 3D positions, a length, an angle, and a torsion.
     *
     * @param x1 A reference to a gemmi::Vec3 object representing the first position.
     * @param x2 A reference to a gemmi::Vec3 object representing the second position.
     * @param x3 A reference to a gemmi::Vec3 object representing the third position.
     * @param length A constant reference to a double representing the length of the projected point.
     * @param angle A constant reference to a double representing the angle of projection.
     * @param torsion A constant reference to a double representing the torsion angle of projection.
     * @return A gemmi::Vec3 object representing the calculated projected point.
     */
    gemmi::Vec3 calculate_projected_point(gemmi::Vec3 &x1, gemmi::Vec3 &x2, gemmi::Vec3 &x3, const double &length,
                                          const double &angle, const double &torsion);


    /**
     * Calculates the superposition of two sets of 3D positions.
     *
     * @param reference A vector of gemmi::Position objects representing the reference positions.
     * @param target A vector of gemmi::Position objects representing the target positions.
     * @return A gemmi::Transform object representing the transformation matrix and translation vector that superimposes the reference positions onto the target positions.
     * @throws std::runtime_error if the size of the reference and target vectors are different.
     */
    gemmi::Transform calculate_superposition(std::vector<gemmi::Position>& reference, std::vector<gemmi::Position>& target);

}

#endif //SAILS_SAILS_VECTOR_H
