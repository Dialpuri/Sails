//
// Created by Jordan Dialpuri on 12/02/2024.
//

#include "sails-model.h"

#include <clipper/core/map_interp.h>


clipper::RTop_orth Sails::LinkageSet::calculate(const clipper::MMonomer& donor_residue,
                                                const clipper::MMonomer& acceptor_residue,
                                                const DonorSet& donor, const AcceptorSet& acceptor,
                                                const LinkageParams& params) {
    const int i1 = donor_residue.lookup(donor.atom1,  clipper::MM::UNIQUE);
    const int i2 = donor_residue.lookup(donor.atom2, clipper::MM::UNIQUE);
    const int i3 = donor_residue.lookup(donor.atom3, clipper::MM::UNIQUE);

    if (i1 < 0 || i2 < 0 || i3 < 0) {
        std::cout << "Found AA without proper atom names, skipping " << i1 << " " << i2 << " " << i3  << std::endl;
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
        std::cout << "Found sugar without proper atom names, aborting" << donor.atom1 << " " << donor.atom2 << " " << donor.atom3  << std::endl;
        exit(-1);
    }
    const clipper::Coord_orth dc1 = acceptor_residue[di1].coord_orth();
    const clipper::Coord_orth dc2 = acceptor_residue[di2].coord_orth();
    const clipper::Coord_orth dc3 = acceptor_residue[di3].coord_orth();

    const clipper::RTop_orth rtop = {{dc1, dc2, dc3}, {position1, position2, position3}};
    return rtop;
}

// Sails::AminoAcid::AminoAcid(clipper::MMonomer& amino_acid,
//                             const std::string& atom_1,
//                             const std::string& atom_2,
//                             const std::string& atom_3,
//                             float phi_angle,
//                             float psi_angle,
//                             float phi_torsion,
//                             float psi_torsion) {
//
//
//     const int i1 = amino_acid.lookup(atom_1,  clipper::MM::UNIQUE);
//     const int i2 = amino_acid.lookup(atom_2, clipper::MM::UNIQUE);
//     const int i3 = amino_acid.lookup(atom_3, clipper::MM::UNIQUE);
//
//     if (i1 < 0 || i2 < 0 || i3 < 0) {
//         std::cout << "Found AA without proper atom names, skipping " << i1 << " " << i2 << " " << i3  << std::endl;
//         return;
//     }
//
//     clipper::Coord_orth c1 = amino_acid[i1].coord_orth();
//     clipper::Coord_orth c2 = amino_acid[i2].coord_orth();
//     clipper::Coord_orth c3 = amino_acid[i3].coord_orth();
//
//     position1 = clipper::Coord_orth(c1,
//                                     c2,
//                                     c3,
//                                     bond_length,
//                                     clipper::Util::d2rad(psi_angle),
//                                     clipper::Util::d2rad(psi_torsion));
//
//     position2 = clipper::Coord_orth(c2,
//                                    c3,
//                                    position1,
//                                    bond_length,
//                                    clipper::Util::d2rad(phi_angle),
//                                    clipper::Util::d2rad(phi_torsion));
//
//     position3 = clipper::Coord_orth(c3,
//                                 position1,
//                                 position2,
//                                 bond_length,
//                                 clipper::Util::d2rad(120),
//                                 clipper::Util::d2rad(180));
//
// }
//
//
// Sails::Sugar::Sugar(clipper::MMonomer& sugar,
//                      const std::string& atom_1,
//                      const std::string& atom_2,
//                      const std::string& atom_3,
//                      clipper::Coord_orth target_1,
//                      clipper::Coord_orth target_2,
//                      clipper::Coord_orth target_3) {
//
//     const int i1 = sugar.lookup(atom_1,  clipper::MM::UNIQUE);
//     const int i2 = sugar.lookup(atom_2, clipper::MM::UNIQUE);
//     const int i3 = sugar.lookup(atom_3, clipper::MM::UNIQUE);
//
//     if (i1 < 0 || i2 < 0 || i3 < 0) {
//         std::cout << "Found sugar without proper atom names, aborting" << std::endl;
//         exit(-1);
//     }
//     clipper::Coord_orth c1 = sugar[i1].coord_orth();
//     clipper::Coord_orth c2 = sugar[i2].coord_orth();
//     clipper::Coord_orth c3 = sugar[i3].coord_orth();
//
//     clipper::RTop_orth rtop = {{c1,c2, c3}, {target_1, target_2, target_3}};
//     sugar.transform(rtop);
//     positioned_monomer = sugar;
// }
//
// float Sails::Sugar::score_sugar(const clipper::Xmap<float>&xmap) {
//
//     if (positioned_monomer.size() == 0) {
//         std::cout << "No atoms in the positioned monomer" << std::endl;
//         return -10;
//     }
//     float score_sum = 0;
//     for (int a = 0; a < positioned_monomer.size(); a++) {
//         clipper::Coord_orth coordinate = positioned_monomer[a].coord_orth();
//         clipper::Coord_frac coordinate_frac = coordinate.coord_frac(xmap.cell());
//         score_sum += xmap.interp<clipper::Interp_cubic>(coordinate_frac);
//     }
//     return score_sum/positioned_monomer.size();
// }

