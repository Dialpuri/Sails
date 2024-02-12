//
// Created by Jordan Dialpuri on 12/02/2024.
//

#ifndef SAILS_MODEL_H
#define SAILS_MODEL_H
#include <clipper/clipper-minimol.h>


class SailsAminoAcid {
public:
    explicit SailsAminoAcid(clipper::MMonomer& aa,
        const std::string& atom_1,
        const std::string& atom_2,
        const std::string& atom_3,
        float angle,
        float torsion);

    virtual ~SailsAminoAcid() = default;

protected:
    clipper::Coord_orth position1;
    clipper::Coord_orth position2;
};

class ASN : public SailsAminoAcid{
public:
    explicit ASN(clipper::MMonomer& aa): SailsAminoAcid(aa, "CB", "CG", "ND2", 110, 180) {};

    clipper::Coord_orth get_position1() {
        return position1;
    }

    clipper::Coord_orth get_position2() {
        return position2;
    }
};

class SailsSugar {
public:
    explicit SailsSugar(clipper::MMonomer& sugar,
        const std::string& atom_1,
        const std::string& atom_2,
        clipper::Coord_orth target_1,
        clipper::Coord_orth target_2);

protected:
    clipper::MMonomer positioned_monomer;
};

class NAG: public SailsSugar {
public:
    explicit NAG(clipper::MMonomer& sugar, clipper::Coord_orth& pos1, clipper::Coord_orth& pos2):
    SailsSugar(sugar, "C1", "O5", pos1, pos2)
{}

    clipper::MMonomer get_positioned_monomer() {
        return positioned_monomer;
    }
};

#endif //SAILS_MODEL_H
