//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_MODEL_H
#define SAILS_SAILS_MODEL_H

#include <map>

#include "gemmi/neighbor.hpp"

namespace Sails {
    /**
     * @class AtomSet
     * @brief A class representing a set of atoms.
     *
     * @constructor AtomSet
     * @param atom_1 Name of atom 1 e.g. C1.
     * @param atom_2 Name of atom 2 e.g. C2.
     * @param atom_3 Name of atom 3 e.g O2.
     * @param identifier Identifier of donor atom e.g. 2
     */
    struct AtomSet {
        AtomSet() = default;

        AtomSet(std::string atom_1, std::string atom_2, std::string atom_3, int identifier)
            : atom1(std::move(atom_1)), atom2(std::move(atom_2)), atom3(std::move(atom_3)), identifier(identifier) {
        }

        [[nodiscard]] std::vector<std::string> get_atom_list() const {
            return {atom1, atom2, atom3};
        }

        std::string atom1;
        std::string atom2;
        std::string atom3;
        int identifier{};
    };


    /**
     * @struct ResidueData
     * @brief A struct representing the data of a residue.
     *
     * This struct stores the acceptor and donor maps, acceptor and donor sets, and SNFG shape
     * and color for a residue.
     *
     * @constructor ResidueData
     * @param acceptors A vector of acceptor atom sets.
     * @param donors A vector of donor atom sets.
     * @param snfg_shape The shape of the SNFG symbol e.g. square.
     * @param snfg_colour The colour of the SNFG symbol e.g. blue.
     */
    struct ResidueData {
        ResidueData() = default;

        ResidueData(const std::vector<AtomSet> &acceptors, const std::vector<AtomSet> &donors, std::string snfg_shape,
                    std::string snfg_colour, std::vector<int> &preferred_depths) : acceptors(acceptors), donors(donors),
            snfg_shape(std::move(snfg_shape)),
            snfg_colour(std::move(snfg_colour)),
            preferred_depths(preferred_depths) {
            for (const auto &acceptor: acceptors) {
                acceptor_map[acceptor.identifier] = acceptor.get_atom_list();
            }

            for (const auto &donor: donors) {
                donor_map[donor.identifier] = donor.get_atom_list();
            }
        }

        std::map<int, std::vector<std::string> > acceptor_map;
        std::map<int, std::vector<std::string> > donor_map;
        std::vector<AtomSet> acceptors;
        std::vector<AtomSet> donors;
        std::string snfg_shape;
        std::string snfg_colour;
        std::vector<int> preferred_depths;
    };

    typedef std::map<std::string, ResidueData> ResidueDatabase;

    /**
     *@struct AngleSet
     *@brief Structure representing alpha, beta and gamma angles
     *
     * @constructor AngleSet
     * @param psi The alpha angle.
     * @param phi The beta angle.
     * @param omega The gamma angle.
     */
    struct AngleSet {
        AngleSet(const double alpha, const double beta, const double gamma) : alpha(alpha), beta(beta), gamma(gamma) {
        }

        std::vector<double> get_in_order() {
            return {alpha, beta, gamma};
        }

        double alpha;
        double beta;
        double gamma;
    };

    /**
     *@struct TorsionAngle
     *@brief Structure representing a single torsion angle with a mean and standard deviation
    *
     * @constructor TorsionAngle
     * @param mean The mean angle.
     * @param stddev The standard deviation of the angle.
     */
    struct TorsionAngle {
        TorsionAngle(const double mean, const double stddev) : mean(mean), stddev(stddev) {
        }

        double mean;
        double stddev;
    };

    /**
     *@struct TorsionSet
     *@brief Structure representing a set of three torsion angles important for representing a glycosidic linkage
     *
     * @constructor TorsionSet
     * @param psi The psi angle.
     * @param phi The phi angle.
     * @param omega The omega angle.
     */
    struct TorsionSet {
        TorsionSet(const TorsionAngle &psi, const TorsionAngle &phi, const TorsionAngle &omega) : psi(psi), phi(phi),
            omega(omega) {
        }

        std::vector<double> get_means_in_order() {
            return {psi.mean, phi.mean, omega.mean};
        }

        TorsionAngle psi;
        TorsionAngle phi;
        TorsionAngle omega;
    };


    /**
     * @struct LinkageData
     * @brief A struct representing the data for a linkage between two atoms.
     *
     * @details The LinkageData struct holds information about a linkage between two atoms, including the donor and acceptor names,
     * the donor and acceptor numbers, and the angle and torsion sets associated with the linkage.
     *
     * @constructor LinkageData
     * @param donor The name of the donor atom.
     * @param acceptor The name of the acceptor atom.
     * @param donor_number The number of the donor atom.
     * @param acceptor_number The number of the acceptor atom.
     * @param angle_set The AngleSet associated with the linkage.
     * @param torsion_set The TorsionSet associated with the linkage.
     */
    struct LinkageData {
        LinkageData(std::string donor, std::string acceptor, const int donor_number, const int acceptor_number,
                    const double length,
                    const AngleSet &angle_set, const TorsionSet &torsion_set) : donor(std::move(donor)),
            acceptor(std::move(acceptor)),
            donor_number(donor_number),
            acceptor_number(acceptor_number), length(length),
            angles(angle_set),
            torsions(torsion_set) {
        };

        std::string donor;
        std::string acceptor;
        int donor_number;
        int acceptor_number;
        double length;
        AngleSet angles;
        TorsionSet torsions;
    };

    typedef std::map<std::string, std::vector<LinkageData> > LinkageDatabase;

    struct Glycosite {
        Glycosite() = default;

        explicit Glycosite(const gemmi::NeighborSearch::Mark &mark) {
            model_idx = 0;
            chain_idx = mark.chain_idx;
            residue_idx = mark.residue_idx;
            atom_idx = mark.atom_idx;
        }

        Glycosite(const int model_idx, const int chain_idx, const int residue_idx) : model_idx(model_idx),
            chain_idx(chain_idx),
            residue_idx(residue_idx) {
        }

        int model_idx{};
        int chain_idx{};
        int residue_idx{};
        int atom_idx = -1; // optional
    };

    typedef std::vector<Glycosite> Glycosites;

    inline std::optional<Glycosite> find_site(const gemmi::Structure &structure,
                                              const std::string &chain_name, const std::string &residue_name,
                                              int seqId) {
        for (int model_idx = 0; model_idx < structure.models.size(); ++model_idx) {
            auto &model = structure.models[model_idx];
            for (int chain_idx = 0; chain_idx < model.chains.size(); ++chain_idx) {
                auto &chain = model.chains[chain_idx];
                if (chain.name != chain_name) continue;
                for (int residue_idx = 0; residue_idx < chain.residues.size(); ++residue_idx) {
                    if (auto &residue = chain.residues[residue_idx]; residue.name == residue_name && residue.seqid.num.value == seqId) {
                        return Glycosite(model_idx, chain_idx, residue_idx);
                    }
                }
            }
        }
        return std::nullopt;
    }
}
#endif //SAILS_SAILS_MODEL_H
