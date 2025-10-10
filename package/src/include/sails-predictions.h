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
        explicit Predictions(gemmi::Grid<>* glycan_map, LinkageDatabase& linkage_database, ResidueDatabase& residue_database): m_residue_database(residue_database) {
            protein_donors = find_protein_donors(linkage_database);
            m_glycan_map = glycan_map;
        };

        explicit Predictions(gemmi::Grid<>* glycan_map, gemmi::Grid<>* protein_map, LinkageDatabase& linkage_database, ResidueDatabase& residue_database):  m_residue_database(residue_database) {
            protein_donors = find_protein_donors(linkage_database);
            m_glycan_map = glycan_map;
            m_protein_map = protein_map;
        };

        Glycosites find_potential_sites(gemmi::Structure &structure);

    private:
        std::optional<gemmi::NeighborSearch> create_neighbour_search(gemmi::Grid<> *grid, float threshold,
                                                                     const gemmi::UnitCell &unit_cell);

        Glycosites find_potential_sites_using_glycan(gemmi::Structure &structure);

        Glycosites find_potential_sites_using_protein_glycan(gemmi::Structure &structure);


        gemmi::Grid<>* m_glycan_map = nullptr;
        gemmi::Grid<>* m_protein_map = nullptr;
        std::set<std::string> protein_donors;
        Sails::ResidueDatabase& m_residue_database;
    };

}


#endif //SAILS_PREDICTIONS_H
