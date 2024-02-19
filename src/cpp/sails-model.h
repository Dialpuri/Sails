//
// Created by Jordan Dialpuri on 12/02/2024.
//
#pragma once
#ifndef SAILS_MODEL_H
#define SAILS_MODEL_H
#include <map>
#include <utility>
#include <clipper/clipper-minimol.h>
#include <clipper/core/xmap.h>

namespace Sails {
    struct LinkageParams {
        LinkageParams() = default;
        LinkageParams(
            const float phi_torsion_angle, const float psi_torsion_angle, const float theta_torsion_angle,
            const float phi_atoms_angle, const float psi_atoms_angle, const float theta_atoms_angle, const float bond_len
        ): phi_torsion(phi_torsion_angle), psi_torsion(psi_torsion_angle), theta_torsion(theta_torsion_angle),
           phi_angle(phi_atoms_angle), psi_angle(psi_atoms_angle), theta_angle(theta_atoms_angle), bond_length(bond_len) {
        }

        float phi_torsion;
        float psi_torsion;
        float theta_torsion;
        float phi_angle;
        float psi_angle;
        float theta_angle;
        float bond_length;

        std::string format() const {
            return " Phi_t = " + std::to_string(phi_torsion) +
                    " Psi_t = " + std::to_string(psi_torsion) +
                    " The_t = " + std::to_string(theta_torsion) +
                    " Phi_a = " + std::to_string(phi_angle) +
                    " Psi_a = " + std::to_string(psi_angle) +
                    " The_a = " + std::to_string(theta_angle);
        }
    };

    struct DonorSet {
        DonorSet() = default;
        DonorSet(const std::string& atom_1, const std::string& atom_2, const std::string& atom_3): atom1(atom_1),
                                                                           atom2(atom_2),
                                                                           atom3(atom_3) {
        }

        std::string atom1;
        std::string atom2;
        std::string atom3;
    };

    struct AcceptorSet {
        AcceptorSet() = default;
        AcceptorSet(const std::string& atom_1, const std::string& atom_2, const std::string& atom_3): atom1(atom_1),
                                                                            atom2(atom_2),
                                                                            atom3(atom_3) {
        }

        std::string atom1;
        std::string atom2;
        std::string atom3;
    };

    struct LinkageSet {
        LinkageSet() = default;
        static clipper::RTop_orth calculate(const clipper::MMonomer&donor_residue,
                                            const clipper::MMonomer&acceptor_residue,
                                            const DonorSet&donor, const AcceptorSet&acceptor,
                                            const LinkageParams&params);
        LinkageParams params;
        virtual clipper::RTop_orth update_transformation(LinkageParams& parameters) {return {};};
        virtual clipper::MMonomer get_donor_monomer() { return {};};
        virtual LinkageParams get_parameters() const { std::cout << "base" << std::endl; return {};}
        virtual void set_parameters(LinkageParams& parameters) { params = parameters;}
    };


    struct NAG  {

        static std::map<int, DonorSet> donors() {
            return {{4, {"C5", "C4", "O4"}}};
        };

        static std::map<int, AcceptorSet> acceptors() {
            return {{1, {"C1", "O5", "C5"}}};
        };
    };

    struct ASN  {
        static std::map<int, DonorSet> donors() {
            return {{1, DonorSet("CB", "CG", "ND2")}};
        };

        static std::map<int, AcceptorSet> acceptors() {
            return {};
        };
    };

    struct ASN_NAG : LinkageSet {
        DonorSet donor = ASN::donors()[1];
        AcceptorSet acceptor = NAG::acceptors()[1];
        LinkageParams params = LinkageParams(-108, 180, 180,
                                            105, 120, 105, 1.45);

        ASN_NAG(const clipper::MMonomer& donor_residue, const clipper::MMonomer& acceptor_residue) {
            m_transformation = calculate(donor_residue, acceptor_residue, donor, acceptor, params);
            m_donor_residue = donor_residue;
            m_acceptor_residue = acceptor_residue;
        }

        clipper::RTop_orth update_transformation(LinkageParams& parameters) override {
            return calculate(m_donor_residue, m_acceptor_residue, donor, acceptor, parameters);
        }

        clipper::RTop_orth get_transformation() const {return m_transformation;}

        clipper::MMonomer get_donor_monomer() override {return m_donor_residue;}

        LinkageParams get_parameters() const override  {return params;}

        void set_parameters(LinkageParams& parameters) override  { params = parameters;}


    private:
        clipper::RTop_orth m_transformation;
        clipper::MMonomer m_donor_residue; 
        clipper::MMonomer m_acceptor_residue;
    
    };

    struct NAG_NAG : LinkageSet {
        DonorSet donor = NAG::donors()[4];
        AcceptorSet acceptor = NAG::acceptors()[1];
        LinkageParams params = LinkageParams(-80, -127, 180,
                                             105, 120, 105, 1.45);

        NAG_NAG(const clipper::MMonomer& donor_residue, const clipper::MMonomer& acceptor_residue) {
            m_transformation = calculate(donor_residue, acceptor_residue, donor, acceptor, params);
            m_donor_residue = donor_residue;
            m_acceptor_residue = acceptor_residue;
        }

        clipper::RTop_orth update_transformation(LinkageParams& parameters) override {
            return calculate(m_donor_residue, m_acceptor_residue, donor, acceptor, parameters);
        }

        [[nodiscard]] clipper::RTop_orth get_transformation() const {return m_transformation;}

        [[nodiscard]] clipper::MMonomer get_donor_monomer() override {return m_donor_residue;}

        [[nodiscard]] LinkageParams get_parameters() const override  {return params;}

        void set_parameters(LinkageParams& parameters) override  { params = parameters;}

    private:
        clipper::RTop_orth m_transformation;
        clipper::MMonomer m_donor_residue;
        clipper::MMonomer m_acceptor_residue;

    };

    enum ResidueType {
        ASN, NAG, BMA, MAN
    };

    static std::map<std::string, ResidueType> ResidueMap = {
            {"ASN", ResidueType::ASN},
            {"NAG", ResidueType::NAG},
            {"BMA", ResidueType::BMA},
            {"MAN", ResidueType::MAN}
    };

}
#endif //SAILS_MODEL_H
