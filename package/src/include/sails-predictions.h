//
// Created by Jordan Dialpuri on 07/10/2025.
//

#ifndef SAILS_PREDICTIONS_H
#define SAILS_PREDICTIONS_H

#include <gemmi/grid.hpp>
#include <gemmi/model.hpp>
#include "sails-model.h"
#include "sails-utils.h"

namespace Sails {

    class Predictions {
    public:
        explicit Predictions(gemmi::Grid<>& glycan_map, LinkageDatabase& linkage_database, ResidueDatabase& residue_database): m_glycan_map(glycan_map), m_residue_database(residue_database) {
            protein_donors = find_protein_donors(linkage_database);
        };

        gemmi::NeighborSearch create_neighbour_search(float threshold, gemmi::UnitCell unit_cell);

        Glycosites find_potential_sites(gemmi::Structure &structure);

    private:
        gemmi::Grid<>& m_glycan_map;
        std::set<std::string> protein_donors;
        Sails::ResidueDatabase& m_residue_database;
    };

}


#endif //SAILS_PREDICTIONS_H
