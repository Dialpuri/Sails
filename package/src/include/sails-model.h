//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_MODEL_H
#define SAILS_SAILS_MODEL_H

#include <iostream>
#include <map>

namespace Sails {
    struct LinkageParams {
        LinkageParams() = default;

        LinkageParams(
                const float phi_torsion_angle, const float psi_torsion_angle, const float theta_torsion_angle,
                const float phi_atoms_angle, const float psi_atoms_angle, const float theta_atoms_angle,
                const float bond_length
        ) : phi_torsion(phi_torsion_angle), psi_torsion(psi_torsion_angle), theta_torsion(theta_torsion_angle),
            phi_angle(phi_atoms_angle), psi_angle(psi_atoms_angle), theta_angle(theta_atoms_angle),
            bond_length(bond_length) {
        }

        float phi_torsion;
        float psi_torsion;
        float theta_torsion;
        float phi_angle;
        float psi_angle;
        float theta_angle;
        float bond_length;

        inline std::string format() const {
            return " Phi_t = " + std::to_string(phi_torsion) +
                   " Psi_t = " + std::to_string(psi_torsion) +
                   " The_t = " + std::to_string(theta_torsion) +
                   " Phi_a = " + std::to_string(phi_angle) +
                   " Psi_a = " + std::to_string(psi_angle) +
                   " The_a = " + std::to_string(theta_angle);
        }
    };

    struct AtomSet {
        AtomSet() = default;

        AtomSet(const std::string &atom_1, const std::string &atom_2, const std::string &atom_3, int identifier) :
                atom1(atom_1),
                atom2(atom_2),
                atom3(atom_3),
                identifier(identifier) {
        }

        std::string atom1;
        std::string atom2;
        std::string atom3;
        int identifier;
    };


    struct ResidueData {

        ResidueData(const std::vector<AtomSet> &acceptors, const std::vector<AtomSet> &donors) :
                acceptors(acceptors),
                donors(donors) {}

        std::vector<AtomSet> acceptors;
        std::vector<AtomSet> donors;
    };

    typedef std::map<std::string, ResidueData> ResidueDatabase;


//    struct LinkageSet {
//        LinkageSet() = default;
//
//        static clipper::RTop_orth calculate(const clipper::MMonomer &donor_residue,
//                                            const clipper::MMonomer &acceptor_residue,
//                                            const AtomSet &donor, const AtomSet &acceptor,
//                                            const LinkageParams &params);
//
//        LinkageParams params;
//
//        virtual clipper::RTop_orth update_transformation(LinkageParams &parameters) { return {}; };
//
//        virtual clipper::MMonomer get_donor_monomer() { return {}; };
//
//        virtual LinkageParams get_parameters() const {
//            std::cout << "base" << std::endl;
//            return {};
//        }
//
//        virtual void set_parameters(LinkageParams &parameters) { params = parameters; }
//    };



}


#endif //SAILS_SAILS_MODEL_H
