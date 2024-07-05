//
// Created by Jordan Dialpuri on 05/07/2024.
//

#include "../include/sails-vector.h"

// Translated from clipper Coord_orth
//const Coord_orth xa( (x3-x2).unit() );
//const Coord_orth xc( Vec3<>::cross( x2-x1, xa ).unit() );
//const Coord_orth xb( Vec3<>::cross( xa, xc ) );
//const ftype wa = -length*cos( angle );
//const ftype wb = -length*sin( angle )*cos( -torsion );
//const ftype wc = -length*sin( angle )*sin( -torsion );
//(*this) = x3 + wa*xa + wb*xb + wc*xc;
//}

gemmi::Vec3 Sails::calculate_projected_point(gemmi::Vec3 &x1, gemmi::Vec3 &x2, gemmi::Vec3 &x3, const double &length,
                                      const double &angle, const double &torsion) {

    const gemmi::Vec3 xa = {x3-x2.normalized()};
    const gemmi::Vec3 xc = {(x2-x1).cross(xa).normalized()};
    const gemmi::Vec3 xb = {xa.cross(xc)};

    const double wa = -length * cos(angle);
    const double wb = -length * sin(angle) * cos(-torsion);
    const double wc = -length * sin(angle) * cos(-torsion);

    return x3 + wa*xa + wb*xb + wc*xc;
}