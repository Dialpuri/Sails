//
// Created by Jordan Dialpuri on 12/02/2024.
//

#include "sails-model.h"



clipper::RTop_orth Sails::LinkageSet::calculate(const clipper::MMonomer& donor_residue,
                                                const clipper::MMonomer& acceptor_residue,
                                                const DonorSet& donor, const AcceptorSet& acceptor,
                                                const LinkageParams& params) {
    const int i1 = donor_residue.lookup(donor.atom1,  clipper::MM::UNIQUE);
    const int i2 = donor_residue.lookup(donor.atom2, clipper::MM::UNIQUE);
    const int i3 = donor_residue.lookup(donor.atom3, clipper::MM::UNIQUE);

    if (i1 < 0 || i2 < 0 || i3 < 0) {
//        std::cout << "Found AA without proper atom names, skipping " << donor.atom1 << " " << donor.atom2 << " " << donor.atom3  << std::endl;
        return {};
    }

    const clipper::Coord_orth c1 = donor_residue[i1].coord_orth();
    const clipper::Coord_orth c2 = donor_residue[i2].coord_orth();
    const clipper::Coord_orth c3 = donor_residue[i3].coord_orth();

    const clipper::Coord_orth position1 = clipper::Coord_orth(c1,
                                    c2,
                                    c3,
                                    params.bond_length,
                                    clipper::Util::d2rad(params.psi_angle),
                                    clipper::Util::d2rad(params.psi_torsion));

    const clipper::Coord_orth position2 = clipper::Coord_orth(c2,
                                   c3,
                                   position1,
                                   params.bond_length,
                                   clipper::Util::d2rad(params.phi_angle),
                                   clipper::Util::d2rad(params.phi_torsion));

    const clipper::Coord_orth position3 = clipper::Coord_orth(c3,
                                position1,
                                position2,
                                params.bond_length,
                                clipper::Util::d2rad(params.theta_angle),
                                clipper::Util::d2rad(params.theta_torsion));

    const int di1 = acceptor_residue.lookup(acceptor.atom1,  clipper::MM::UNIQUE);
    const int di2 = acceptor_residue.lookup(acceptor.atom2, clipper::MM::UNIQUE);
    const int di3 = acceptor_residue.lookup(acceptor.atom3, clipper::MM::UNIQUE);

    if (di1 < 0 || di2 < 0 || di3 < 0) {
        std::cout << "Found sugar without proper atom names, aborting" << acceptor.atom1 << " " << acceptor.atom2 << " " << acceptor.atom3  << std::endl;
        exit(-1);
    }
    const clipper::Coord_orth dc1 = acceptor_residue[di1].coord_orth();
    const clipper::Coord_orth dc2 = acceptor_residue[di2].coord_orth();
    const clipper::Coord_orth dc3 = acceptor_residue[di3].coord_orth();

    const clipper::RTop_orth rtop = {{dc1, dc2, dc3}, {position1, position2, position3}};
    return rtop;
}
