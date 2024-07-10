//
// Created by Jordan Dialpuri on 05/07/2024.
//

#include "../include/sails-vector.h"


// Translated from clipper Coord_orth
gemmi::Vec3 Sails::calculate_projected_point(gemmi::Vec3 &x1, gemmi::Vec3 &x2, gemmi::Vec3 &x3, const double &length,
                                      const double &angle, const double &torsion) {

    const gemmi::Vec3 xa = {(x3-x2).normalized()};
    const gemmi::Vec3 xc = {(x2-x1).cross(xa).normalized()};
    const gemmi::Vec3 xb = {xa.cross(xc)};

    const double wa = -length * cos(angle);
    const double wb = -length * sin(angle) * cos(-torsion);
    const double wc = -length * sin(angle) * sin(-torsion);

    return x3 + wa*xa + wb*xb + wc*xc;
}

// Translated from clipper RTop_orth
gemmi::Transform Sails::calculate_superposition(std::vector<gemmi::Position> &source, std::vector<gemmi::Position> &target) {

    if ( source.size() != target.size() )
        throw std::runtime_error("RTop_orth: coordinate list size mismatch");

    // get centre of mass
    clipper::Coord_orth s, t, p, m;
    clipper::Coord_orth src_cen(0.0,0.0,0.0);
    clipper::Coord_orth tgt_cen(0.0,0.0,0.0);
    std::vector<clipper::Coord_orth> source_clipper;
    std::vector<clipper::Coord_orth> target_clipper;

    int n = source.size();
    for ( int i = 0; i < n; i++ ) {
        auto _s = clipper::Coord_orth(source[i].x, source[i].y, source[i].z);
        auto _t = clipper::Coord_orth(target[i].x, target[i].y, target[i].z);

        src_cen = src_cen + _s;
        tgt_cen = tgt_cen + _t;
        source_clipper.emplace_back(_s);
        target_clipper.emplace_back(_t);
    }
    src_cen = (1.0/n) * src_cen;
    tgt_cen = (1.0/n) * tgt_cen;

    // prepare cross-sums
    clipper::Matrix<> mat( 4, 4, 0.0 );
    for ( int i = 0; i < n; i++ ) {
     s = source_clipper[i] - src_cen;
     t = target_clipper[i] - tgt_cen;
     p = s + t;
     m = s - t;
     mat(0,0) = mat(0,0) + m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
     mat(1,1) = mat(1,1) + m[0]*m[0] + p[1]*p[1] + p[2]*p[2];
     mat(2,2) = mat(2,2) + p[0]*p[0] + m[1]*m[1] + p[2]*p[2];
     mat(3,3) = mat(3,3) + p[0]*p[0] + p[1]*p[1] + m[2]*m[2];
     mat(0,1) = mat(0,1) + m[2]*p[1] - m[1]*p[2];
     mat(0,2) = mat(0,2) + m[0]*p[2] - m[2]*p[0];
     mat(0,3) = mat(0,3) + m[1]*p[0] - m[0]*p[1];
     mat(1,2) = mat(1,2) + m[0]*m[1] - p[0]*p[1];
     mat(1,3) = mat(1,3) + m[0]*m[2] - p[0]*p[2];
     mat(2,3) = mat(2,3) + m[1]*m[2] - p[1]*p[2];
    }
    mat(1,0) = mat(0,1);
    mat(2,0) = mat(0,2);
    mat(2,1) = mat(1,2);
    mat(3,0) = mat(0,3);
    mat(3,1) = mat(1,3);
    mat(3,2) = mat(2,3);
    // eigenvalue calc
    std::vector<double> v = mat.eigen( true );
    // get result
    clipper::Rotation r( mat(0,0), mat(1,0), mat(2,0), mat(3,0) );
    clipper::Mat33<> rot = r.norm().matrix();
    auto rtop = clipper::RTop_orth( rot, tgt_cen - rot*src_cen );
    auto translation_clipper = rtop.trn();
    auto rotation_clipper = rtop.rot();
    gemmi::Vec3 translation = {translation_clipper[0], translation_clipper[1], translation_clipper[2]};
    gemmi::Mat33 matrix = {
        rotation_clipper(0,0), rotation_clipper(0,1), rotation_clipper(0, 2),
        rotation_clipper(1,0), rotation_clipper(1,1), rotation_clipper(1, 2),
        rotation_clipper(2,0), rotation_clipper(2,1), rotation_clipper(2, 2)
    };
    gemmi::Transform transform;
    transform.mat = matrix;
    transform.vec = translation;

    return transform;
}