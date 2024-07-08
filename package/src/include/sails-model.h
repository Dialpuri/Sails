//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_MODEL_H
#define SAILS_SAILS_MODEL_H

#include <iostream>
#include <map>

#include "gemmi/neighbor.hpp"

namespace Sails {
    struct LinkageParams {
        LinkageParams() = default;

        LinkageParams(const float phi_torsion_angle, const float psi_torsion_angle, const float theta_torsion_angle,
                      const float phi_atoms_angle, const float psi_atoms_angle, const float theta_atoms_angle,
                      const float bond_length) : phi_torsion(phi_torsion_angle), psi_torsion(psi_torsion_angle),
                                                 theta_torsion(theta_torsion_angle), phi_angle(phi_atoms_angle),
                                                 psi_angle(psi_atoms_angle), theta_angle(theta_atoms_angle),
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
            return " Phi_t = " + std::to_string(phi_torsion) + " Psi_t = " + std::to_string(psi_torsion) + " The_t = " +
                   std::to_string(theta_torsion) + " Phi_a = " + std::to_string(phi_angle) + " Psi_a = " +
                   std::to_string(psi_angle) + " The_a = " + std::to_string(theta_angle);
        }
    };

    struct AtomSet {
        AtomSet() = default;

        AtomSet(const std::string &atom_1, const std::string &atom_2, const std::string &atom_3, int identifier)
                : atom1(atom_1), atom2(atom_2), atom3(atom_3), identifier(identifier) {
        }

        std::vector<std::string> get_atom_list() const {
            return {atom1, atom2, atom3};
        }

        std::string atom1;
        std::string atom2;
        std::string atom3;
        int identifier;
    };


    struct ResidueData {
        ResidueData() = default;

        ResidueData(std::vector<AtomSet> &acceptors, std::vector<AtomSet> &donors, const std::string &snfg_shape,
                    const std::string &snfg_colour) : acceptors(acceptors), donors(donors), snfg_shape(snfg_shape),
                                                      snfg_colour(snfg_colour) {
            for (const auto &acceptor: acceptors) {
                acceptor_map[acceptor.identifier] = acceptor.get_atom_list();
            }

            for (const auto &donor: donors) {
                donor_map[donor.identifier] = donor.get_atom_list();
            }
        }

        std::map<int, std::vector<std::string>> acceptor_map;
        std::map<int, std::vector<std::string>> donor_map;
        std::vector<AtomSet> acceptors;
        std::vector<AtomSet> donors;
        std::string snfg_shape;
        std::string snfg_colour;
    };

    typedef std::map<std::string, ResidueData> ResidueDatabase;

    struct AngleSet {
        AngleSet(double alpha, double beta, double gamma) : alpha(alpha), beta(beta), gamma(gamma) {}

        std::vector<double> get_in_order() {
            return {alpha, beta, gamma};
        }

        double alpha;
        double beta;
        double gamma;
    };

    struct TorsionAngle {
        TorsionAngle(double mean, double stddev) : mean(mean), stddev(stddev) {}

        double mean;
        double stddev;
    };

    struct TorsionSet {
        TorsionSet(TorsionAngle &psi, TorsionAngle &phi, TorsionAngle &omega) : psi(psi), phi(phi), omega(omega) {}

        std::vector<double> get_means_in_order() {
            return {psi.mean, phi.mean, omega.mean};
        }

        TorsionAngle psi;
        TorsionAngle phi;
        TorsionAngle omega;
    };

    struct LinkageData {
        LinkageData(const std::string &donor, const std::string &acceptor, int donor_number, int acceptor_number,
                    AngleSet &angle_set, TorsionSet &torsion_set) : donor(donor), acceptor(acceptor),
                                                                    donor_number(donor_number),
                                                                    acceptor_number(acceptor_number), angles(angle_set),
                                                                    torsions(torsion_set) {};

        std::string donor;
        std::string acceptor;
        int donor_number;
        int acceptor_number;
        AngleSet angles;
        TorsionSet torsions;

    };

    typedef std::map<std::string, std::vector<LinkageData>> LinkageDatabase;

    struct Glycosite {
        Glycosite() = default;

        explicit Glycosite(gemmi::NeighborSearch::Mark &mark) {
            model_idx = 0;
            chain_idx = mark.chain_idx;
            residue_idx = mark.residue_idx;
            atom_idx = mark.atom_idx;
        }

        Glycosite(int model_idx, int chain_idx, int residue_idx) : model_idx(model_idx), chain_idx(chain_idx),
                                                                   residue_idx(residue_idx) {}

        int model_idx;
        int chain_idx;
        int residue_idx;
        int atom_idx = -1; // optional
    };

    typedef std::vector<Glycosite> Glycosites;

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
