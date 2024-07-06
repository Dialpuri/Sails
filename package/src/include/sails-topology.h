//
// Created by Jordan Dialpuri on 06/07/2024.
//

#ifndef SAILS_SAILS_TOPOLOGY_H
#define SAILS_SAILS_TOPOLOGY_H

#include "sails-glycan.h"
#include "sails-model.h"

#include "gemmi/model.hpp"
#include "gemmi/neighbor.hpp"

namespace Sails {

    std::optional<Sails::Glycan> find_glycan_topology(gemmi::Structure& structure, const Sails::Glycosite& glycosite);


}

#endif //SAILS_SAILS_TOPOLOGY_H
