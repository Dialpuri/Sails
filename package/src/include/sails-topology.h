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

    struct Topology {
        Topology(gemmi::Structure& structure, Sails::ResidueDatabase& database);

        Glycan find_glycan_topology(Sails::Glycosite &glycosite);


    private:
        gemmi::Structure m_structure;
        ResidueDatabase m_database;
        gemmi::NeighborSearch m_neighbor_search;

        void find_residue_near_donor(Sails::Glycosite &glycosite,  Sails::Glycan &glycan, std::queue<Glycosite> &queue);

    };
}

#endif //SAILS_SAILS_TOPOLOGY_H
