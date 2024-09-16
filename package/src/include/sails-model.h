//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_MODEL_H
#define SAILS_SAILS_MODEL_H

#include <map>
#include <optional>

#include <gemmi/neighbor.hpp>

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
     * @param preferred_depths A vector of preferred depths for the sugar, e.g. NAG prefers 1 and 2.
     * @param anomer The anomeric designation (alpha or beta)
     * @param special Whether the residue needs special treatment in the SNFG
     */
    struct ResidueData {
        ResidueData() = default;

        ResidueData(const std::vector<AtomSet> &acceptors, const std::vector<AtomSet> &donors, std::string& snfg_shape,
                    std::string& snfg_colour, std::vector<int> &preferred_depths, std::string& anomer,
                    std::string& wurcs
                    ) : acceptors(acceptors), donors(donors),
            snfg_shape(std::move(snfg_shape)),
            snfg_colour(std::move(snfg_colour)),
            preferred_depths(preferred_depths),
            anomer(anomer),
            special(special) {

            if (!wurcs.empty()) {wurcs_code = wurcs;}

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
        std::string anomer;
        bool special;
        std::optional<std::string> wurcs_code = std::nullopt;
    };

    typedef std::map<std::string, ResidueData> ResidueDatabase;

    /**
        *@struct Angle
        *@brief Structure representing a single angle with a mean and standard deviation
       *
        * @constructor Angle
        * @param mean The mean angle.
        * @param stddev The standard deviation of the angle.
        */
    struct Angle {
        Angle(const double mean, const double stddev) : mean(mean), stddev(stddev) {
        }

        double mean;
        double stddev;
    };

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
        AngleSet(const Angle &alpha, const Angle &beta, const Angle &gamma) : alpha(alpha), beta(beta), gamma(gamma) {
        }

        std::vector<double> get_means_in_order() {
            return {alpha.mean, beta.mean, gamma.mean};
        }

        std::vector<double> get_stddev_in_order() {
            return {alpha.stddev, beta.stddev, gamma.stddev};
        }

        Angle alpha;
        Angle beta;
        Angle gamma;
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
        TorsionSet(const Angle &psi, const Angle &phi, const Angle &omega) : psi(psi), phi(phi),
                                                                             omega(omega) {
        }

        std::vector<double> get_means_in_order() {
            return {psi.mean, phi.mean, omega.mean};
        }

        std::vector<double> get_stddev_in_order() {
            return {psi.stddev, phi.stddev, omega.stddev};
        }

        Angle psi;
        Angle phi;
        Angle omega;
    };


    /**
     * @class Cluster
     * @brief A class representing a cluster of angles and torsions.
     *
     * @constructor Cluster
     * @param angles A set of angles.
     * @param torsions A set of torsions.
     */
    struct Cluster {
        Cluster(const AngleSet &angles, const TorsionSet &torsions)
            : angles(angles),
              torsions(torsions) {
        }

        AngleSet angles;
        TorsionSet torsions;
    };

    typedef std::vector<Cluster> Clusters;
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
                    const std::vector<Cluster> &clusters) : donor(std::move(donor)),
                                                            acceptor(std::move(acceptor)),
                                                            donor_number(donor_number),
                                                            acceptor_number(acceptor_number), length(length),
                                                            clusters(clusters) {
        };

        std::string donor;
        std::string acceptor;
        int donor_number;
        int acceptor_number;
        double length;
        std::vector<Cluster> clusters;
    };

    typedef std::map<std::string, std::vector<LinkageData> > LinkageDatabase;

    /**
     * @class Glycosite
     * @brief A class representing a glycosite.
     *
     * The Glycosite class represents a glycosite, which is defined by its model index, chain index, residue index,
     * and optional atom index. Glycosite can be used for set operations such as comparison and sorting.
     *
     *
     * @fn bool Glycosite::operator<(const Glycosite &other) const
     * @brief Less than operator.
     * @param other The Glycosite object to compare with.
     * @return True if this Glycosite object is less than the other Glycosite object, False otherwise.
     *
     * Compares this Glycosite object with the other Glycosite object based on model_idx, chain_idx,
     * residue_idx, and atom_idx. Returns True if this Glycosite object is less than the other Glycosite
     * object, and False otherwise.
     *
     * @var int Glycosite::model_idx
     * @brief The model index of the glycosite.
     *
     * This variable stores the model index of the glycosite.
     *
     * @var int Glycosite::chain_idx
     * @brief The chain index of the glycosite.
     *
     * This variable stores the chain index of the glycosite.
     *
     * @var int Glycosite::residue_idx
     * @brief The residue index of the glycosite.
     *
     * This variable stores the residue index of the glycosite.
     *
     * @var int Glycosite::atom_idx
     * @brief The atom index of the glycosite.
     *
     * This variable stores the atom index of the glycosite. It is set to -1 by default and is optional.
     */
    struct Glycosite {
        /**
         * @constructor Glycosite
         * @brief Default constructor.
         *
         * Constructs a Glycosite object with default values for model_idx, chain_idx, residue_idx, and atom_idx.
         *
         */
        Glycosite() = default;


        /**
         * @constructor Glycosite
         * @brief Constructor with mark parameter.
         * @param mark The neighbor search mark used to initialize the Glycosite object.
         *
         * Constructs a Glycosite object using the provided neighbor search mark. The model_idx is set to 0,
         * and the chain_idx, residue_idx, and atom_idx are set based on the mark.
         *
         */
        explicit Glycosite(const gemmi::NeighborSearch::Mark &mark) {
            model_idx = 0;
            chain_idx = mark.chain_idx;
            residue_idx = mark.residue_idx;
            atom_idx = mark.atom_idx;
        }

        /**
         * @constructor Glycosite
         * @brief Constructor with model_idx, chain_idx, and residue_idx parameters.
         * @param model_idx The model index of the glycosite.
         * @param chain_idx The chain index of the glycosite.
         * @param residue_idx The residue index of the glycosite.
         *
         * Constructs a Glycosite object using the provided model index, chain index, and residue index.
         * The atom_idx is set to -1 by default.
         */
        Glycosite(const int model_idx, const int chain_idx, const int residue_idx) : model_idx(model_idx),
            chain_idx(chain_idx),
            residue_idx(residue_idx) {
        }

        /**
         * @constructor Glycosite
         * @brief Constructor with model_idx, chain_idx, and residue_idx parameters.
         * @param model_idx The model index of the glycosite.
         * @param chain_idx The chain index of the glycosite.
         * @param residue_idx The residue index of the glycosite.
         * @param atom_idx The atom index of the glycosite.
         *
         * Constructs a Glycosite object using the provided model index, chain index, and residue index.
         */
        Glycosite(const int model_idx, const int chain_idx, const int residue_idx,
                  const int atom_idx) : model_idx(model_idx),
                                        chain_idx(chain_idx),
                                        residue_idx(residue_idx),
                                        atom_idx(atom_idx) {
        }

        /**
         * @brief Overloaded less than operator for the Glycosite struct.
         *
         * This operator compares two Glycosite structs based on their model_idx, chain_idx, residue_idx, and atom_idx
         * members. It returns true if this Glycosite is less than the other Glycosite, and false otherwise.
         *
         * @param other The Glycosite struct to compare.
         * @return true if this Glycosite is less than the other Glycosite, false otherwise.
         */
        bool operator<(const Glycosite &other) const {
            return std::tie(model_idx, chain_idx, residue_idx, atom_idx) < std::tie(
                       other.model_idx, other.chain_idx, other.residue_idx, other.atom_idx);
        }

        bool operator==(const Glycosite &other) const {
            return std::tie(model_idx, chain_idx, residue_idx, atom_idx) == std::tie(
                       other.model_idx, other.chain_idx, other.residue_idx, other.atom_idx);
        }


        int model_idx{};
        int chain_idx{};
        int residue_idx{};
        int atom_idx = -1; // optional
    };

    typedef std::vector<Glycosite> Glycosites;

    /**
     * @brief Finds a glycosite in a given structure.
     *
     * The method searches for a glycosite in the provided structure based on the specified chain name,
     * residue name, and sequence ID. It returns an optional Glycosite object if found, otherwise it
     * returns an empty optional.
     *
     * @param structure The input structure to search for the glycosite.
     * @param chain_name The name of the chain where the glycosite is located.
     * @param residue_name The name of the residue that represents the glycosite.
     * @param seqId The sequence identifier of the residue in the chain.
     * @return An optional Glycosite object representing the found glycosite, or an empty optional
     *         if no glycosite is found.
     */
    inline std::optional<Glycosite> find_site(const gemmi::Structure &structure,
                                              const std::string &chain_name, const std::string &residue_name,
                                              int seqId) {
        for (int model_idx = 0; model_idx < structure.models.size(); ++model_idx) {
            auto &model = structure.models[model_idx];
            for (int chain_idx = 0; chain_idx < model.chains.size(); ++chain_idx) {
                auto &chain = model.chains[chain_idx];
                if (chain.name != chain_name) continue;
                for (int residue_idx = 0; residue_idx < chain.residues.size(); ++residue_idx) {
                    if (auto &residue = chain.residues[residue_idx];
                        residue.name == residue_name && residue.seqid.num.value == seqId) {
                        return Glycosite(model_idx, chain_idx, residue_idx);
                    }
                }
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Finds a glycosite in a given structure.
     *
     * The method searches for a glycosite in the provided structure based on the specified chain name,
     * and sequence ID. It returns an optional Glycosite object if found, otherwise it
     * returns an empty optional.
     *
     * @param structure The input structure to search for the glycosite.
     * @param chain_name The name of the chain where the glycosite is located.
     * @param seqId The sequence identifier of the residue in the chain.
     * @return An optional Glycosite object representing the found glycosite, or an empty optional
     *         if no glycosite is found.
     */
    inline std::optional<Glycosite> find_site(const gemmi::Structure &structure,
                                              const std::string &chain_name,
                                              int seqId) {
        for (int model_idx = 0; model_idx < structure.models.size(); ++model_idx) {
            auto &model = structure.models[model_idx];
            for (int chain_idx = 0; chain_idx < model.chains.size(); ++chain_idx) {
                auto &chain = model.chains[chain_idx];
                if (chain.name != chain_name) continue;
                for (int residue_idx = 0; residue_idx < chain.residues.size(); ++residue_idx) {
                    if (auto &residue = chain.residues[residue_idx];
                        residue.seqid.num.value == seqId) {
                        return Glycosite(model_idx, chain_idx, residue_idx);
                    }
                }
            }
        }
        return std::nullopt;
    }
}
#endif //SAILS_SAILS_MODEL_H
