//
// Created by Jordan Dialpuri on 05/07/2024.
//

#ifndef SAILS_SAILS_VECTOR_H
#define SAILS_SAILS_VECTOR_H

#include "gemmi/math.hpp"

namespace Sails {

    gemmi::Vec3 calculate_projected_point(gemmi::Vec3 &x1, gemmi::Vec3 &x2, gemmi::Vec3 &x3, const double &length,
                                                 const double &angle, const double &torsion);
}

#endif //SAILS_SAILS_VECTOR_H
