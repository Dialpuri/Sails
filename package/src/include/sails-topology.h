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

    std::optional<Sails::Glycan>
    find_glycan_topology(gemmi::Structure &structure, Sails::Glycosite &glycosite, Sails::ResidueDatabase &database);

    void find_residue_near_donor(Sails::Glycosite &glycosite, gemmi::Structure &structure,
                                 gemmi::NeighborSearch &neighbor_search, Sails::Glycan &glycan,
                                 Sails::ResidueDatabase &database, std::queue<Glycosite> &queue);


}

#endif //SAILS_SAILS_TOPOLOGY_H
