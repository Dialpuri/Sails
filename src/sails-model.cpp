//
// Created by Jordan Dialpuri on 12/02/2024.
//

#include "sails-model.h"


SailsAminoAcid::SailsAminoAcid(clipper::MMonomer& aa,
                                const std::string& atom_1,
                                const std::string& atom_2,
                                const std::string& atom_3,
                                float angle,
                                float torsion) {

    std::cout << atom_1 << " " << atom_2 << " " << atom_3 << std::endl;
    const int i1 = aa.lookup(atom_1,  clipper::MM::UNIQUE);
    const int i2 = aa.lookup(atom_2, clipper::MM::UNIQUE);
    const int i3 = aa.lookup(atom_3, clipper::MM::UNIQUE);

    if (i1 < 0 || i2 < 0 || i3 < 0) {
        std::cout << "Found AA without proper atom names, aborting" << std::endl;
        exit(-1);
    }

    clipper::Coord_orth c1 = aa[i1].coord_orth();
    clipper::Coord_orth c2 = aa[i2].coord_orth();
    clipper::Coord_orth c3 = aa[i3].coord_orth();

    position1 = clipper::Coord_orth(c1,
                                    c2,
                                    c3,
                                    1.45,
                                    clipper::Util::d2rad(angle),
                                    clipper::Util::d2rad(torsion));
    position2 = clipper::Coord_orth(c2,
                                   c3,
                                   position1,
                                   1.45,
                                   clipper::Util::d2rad(105),
                                   clipper::Util::d2rad(-86));

}

SailsSugar::SailsSugar(clipper::MMonomer& sugar,
                        const std::string& atom_1,
                        const std::string& atom_2,
                        clipper::Coord_orth target_1,
                        clipper::Coord_orth target_2) {

    const int i1 = sugar.lookup(atom_1,  clipper::MM::UNIQUE);
    const int i2 = sugar.lookup(atom_2, clipper::MM::UNIQUE);
    if (i1 < 0 || i2 < 0) {
        std::cout << "Found sugar without proper atom names, aborting" << std::endl;
        exit(-1);
    }
    clipper::Coord_orth c1 = sugar[i1].coord_orth();
    clipper::Coord_orth c2 = sugar[i2].coord_orth();

    clipper::RTop_orth rtop = {{c1,c2}, {target_1, target_2}};
    sugar.transform(rtop);
    positioned_monomer = sugar;
}

