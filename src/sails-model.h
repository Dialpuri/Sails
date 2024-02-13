//
// Created by Jordan Dialpuri on 12/02/2024.
//

#ifndef SAILS_MODEL_H
#define SAILS_MODEL_H
#include <map>
#include <utility>
#include <clipper/clipper-minimol.h>
#include <clipper/core/xmap.h>

namespace Sails {
    struct LinkageParams {
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
        static clipper::RTop_orth calculate(const clipper::MMonomer&donor_residue,
                                            const clipper::MMonomer&acceptor_residue,
                                            const DonorSet&donor, const AcceptorSet&acceptor,
                                            const LinkageParams&params);
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
            transformation = calculate(donor_residue, acceptor_residue, donor, acceptor, params);
        }

        clipper::RTop_orth get_transformation() const {return transformation;}

    private:
        clipper::RTop_orth transformation;
    };


    enum AminoAcidType {
        ASN
    };

    static std::map<std::string, AminoAcidType> AminoAcidMap = {{"ASN", AminoAcidType::ASN}};
    //
    // /**
    //  * \class AminoAcid
    //  * \brief Represents an amino acid in Sails.
    //  *
    //  * The SailsAminoAcid class contains the information about an amino acid in Sails.
    //  * After construction, it calculates the the two joining sugar atom positions defined by the psi and phi angles of an ASN-pyranose linkage
    //  * https://doi.org/10.1107%2FS2059798323003510
    //  */
    // class AminoAcid {
    // public:
    //     /**
    //      * @brief SailsAminoAcid constructor.
    //      * @param amino_acid The reference to the amino acid MMonomer.
    //      * @param atom_1 The name of the first atom in the phi/psi definition.
    //      * @param atom_2 The name of the second atom in the phi/psi definition.
    //      * @param atom_3 The name of the third atom in the phi/psi definition.
    //      * @param phi_angle The angle value between the terminal aa atom, first atom of sugar and second atom of sugar.
    //      * @param psi_angle The angle value between the penultimate aa atom, terminal aa atom and first atom of sugar.
    //      * @param phi_torsion The torsion value phi torsion angle.
    //      * @param psi_torsion The torsion value psi torsion angle.
    //      */
    //     explicit AminoAcid(clipper::MMonomer&amino_acid,
    //                        const std::string&atom_1,
    //                        const std::string&atom_2,
    //                        const std::string&atom_3,
    //                        float phi_angle,
    //                        float psi_angle,
    //                        float phi_torsion,
    //                        float psi_torsion);
    //
    //     virtual ~AminoAcid() = default;
    //
    // protected:
    //     /**
    //      * @brief First atom of pyranose sugar that would be linked to this amino acid, for ASN-NAG position1 would represent C1.
    //      */
    //     clipper::Coord_orth position1;
    //     /**
    //     * @brief Second atom of pyranose sugar that would be linked to this amino acid, for ASN-NAG position2 would represent O5.
    //     */
    //     clipper::Coord_orth position2;
    //     /**
    //     * @brief Second atom of pyranose sugar that would be linked to this amino acid, for ASN-NAG position2 would represent O5.
    //     */
    //     clipper::Coord_orth position3;
    //
    //
    //     float bond_length = 1.45;
    // };
    //
    // /**
    //  * @class Asparagine
    //  * @brief A class representing the ASN amino acid
    //  * Inherits from SailsAminoAcid class.
    //  */
    // class Asparagine : public AminoAcid {
    // public:
    //     /**
    //      * @brief Constructor for ASN class.
    //      * @param aa The reference to the clipper::MMonomer object.
    //      * @details This constructor initializes a new ASN object by calling the base class constructor SailsAminoAcid.
    //      * It sets the atom_1, atom_2, atom_3, angle, and torsion parameters to "CB", "CG", "ND2", 110, and 180 respectively.
    //      */
    //     explicit Asparagine(clipper::MMonomer&aa): AminoAcid(aa, "CB", "CG", "ND2", 105, 120, -86, 180) {
    //     };
    //
    //     /**
    //      * @brief Get the value of position1.
    //      * @see SailsAminoAcid
    //      * @return The position1 value of type Coord_orth.
    //      */
    //     clipper::Coord_orth get_position1() const {
    //         return position1;
    //     }
    //
    //     /**
    //      * @brief Get the value of position2.
    //      * @see SailsAminoAcid
    //      * @return The position2 coordinate of the Coord_orth.
    //      */
    //     clipper::Coord_orth get_position2() const {
    //         return position2;
    //     }
    //
    //     /**
    //     * @brief Get the value of position3.
    //     * @see SailsAminoAcid
    //     * @return The position2 coordinate of the Coord_orth.
    //     */
    //     clipper::Coord_orth get_position3() const {
    //         return position3;
    //     }
    // };
    //
    // /**
    //  * @class Sugar
    //  * @brief A class representing a sugar molecule in a Sails model.
    //  */
    // class Sugar {
    // public:
    //     /**
    //      * @brief Constructs a Sugar object.
    //      *
    //      * @param sugar A reference to a MMonomer object representing the sugar molecule.
    //      * @param atom_1 A std::string representing the name of the first atom.
    //      * @param atom_2 A std::string representing the name of the second atom.
    //      * @param target_1 A Coord_orth object representing the target coordinate of the first atom.
    //      * @param target_2 A Coord_orth object representing the target coordinate of the second atom.
    //      *
    //      * @return None
    //      */
    //     explicit Sugar(clipper::MMonomer&sugar,
    //                    const std::string&atom_1,
    //                    const std::string&atom_2,
    //                    const std::string&atom_3,
    //                    clipper::Coord_orth target_1,
    //                    clipper::Coord_orth target_2,
    //                    clipper::Coord_orth target_3);
    //
    //     float score_sugar(const clipper::Xmap<float>&xmap);
    //
    // protected:
    //     /**
    //      * @brief The transformed sugar into the position defined my target1 and target2
    //      */
    //     clipper::MMonomer positioned_monomer;
    // };
    //
    // /**
    //  * @class NAG
    //  * @brief Represents a NAG (N-Acetylglucosamine) molecule.
    //  *
    //  * This class inherits from SailsSugar and provides a specialized implementation
    //  * for NAG molecules. It includes methods to get the positioned monomer.
    //  */
    // class NAG : public Sugar {
    // public:
    //     /**
    //      * @brief Constructs a NAG object using a sugar monomer, and two orthogonal coordinates positions
    //      * @param sugar The sugar monomer object
    //      * @param pos1 The first orthogonal coordinate position
    //      * @param pos2 The second orthogonal coordinate position
    //      */
    //     explicit NAG(clipper::MMonomer&sugar, clipper::Coord_orth&pos1, clipper::Coord_orth&pos2,
    //                  clipper::Coord_orth&pos3): Sugar(sugar, "C1", "O5", "C5", pos1, pos2, pos3) {
    //     }
    //
    //     /**
    //      * @brief Retrieves the positioned monomer.
    //      *
    //      * This method returns the positioned monomer object of the SailsSugar class.
    //      *
    //      * @return The positioned monomer object.
    //      */
    //     clipper::MMonomer get_positioned_monomer() {
    //         return positioned_monomer;
    //     }
    // };
}
#endif //SAILS_MODEL_H
