//
// Created by Jordan Dialpuri on 06/07/2024.
//

#ifndef SAILS_SAILS_TOPOLOGY_H
#define SAILS_SAILS_TOPOLOGY_H

#include "sails-glycan.h"
#include "sails-utils.h"

#include "gemmi/model.hpp"
#include "gemmi/neighbor.hpp" // for neighboursearch
#include "gemmi/resinfo.hpp" // for find_tabulated_resiude

namespace Sails {
    /**
     * @class Topology
     * @brief Class for finding glycan topology in a structure
     *
     * The Topology class is responsible for finding the glycan topology in a given structure by using a provided residue database
     * and looking for donor-acceptor atoms within a specified residue.
     *
     * @constructor Topology
     * @param structure The structure to search
     * @param database The residue database which contains information about expected atom donors and acceptors.
     */
    struct Topology {

        Topology(gemmi::Structure* structure, Sails::ResidueDatabase& database);

        /**
         * @brief Find the glycan topology in a given glycosite
         *
         * This method is responsible for finding the glycan topology in a given glycosite by using a provided structure and residue database.
         * It iterates through each residue in the glycosite and checks for close donor-acceptor atom pairs.
         *
         * @param glycosite The glycosite to search in
         * @return The glycan topology found in the glycosite
         */
        Glycan find_glycan_topology(Sails::Glycosite &glycosite);

        void set_structure(gemmi::Structure* structure) {
            m_structure = structure;
            initialise_neighbour_search(structure);
        }

        /**
         * @brief Get the structure object.
         *
         * This method returns the internal structure object.
         *
         * @return The gemmi::Structure object.
         */
        [[nodiscard]] gemmi::Structure get_structure() const {
            return *m_structure;
        }

    private:

        void initialise_neighbour_search(gemmi::Structure* structure);


        /**
         * @brief Find nearby residue near a donor atom.
         *
         * This method finds nearby residues for a given donor atom in a glycosite. It adds the found residues to the glycan and adds the corresponding linkages to the glycan as well.
         *
         * @param glycosite The glycosite containing the donor atom.
         * @param glycan The glycan to add the found residues and linkages to.
         * @param queue The queue to store the found glycosites for further processing.
         *
         * @throws std::runtime_error If the glycosite is not found in the database.
         */
        void find_residue_near_donor(Sails::Glycosite &glycosite,  Sails::Glycan &glycan, std::queue<Glycosite> &queue);

    private:
        gemmi::Structure* m_structure;
        ResidueDatabase m_database;
        gemmi::NeighborSearch m_neighbor_search;

    };
}

#endif //SAILS_SAILS_TOPOLOGY_H
