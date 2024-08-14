//
// Created by Jordan Dialpuri on 13/08/2024.
//

#include "../../include/density/sails-density.h"
#include "../../include/density/sails-em-density.h"

Sails::EMDensity::EMDensity(gemmi::Grid<> &grid) {
    m_grid = grid;

    gemmi::Ccp4<> m;
    m.grid = grid;
    m.update_ccp4_header();
    m.write_ccp4_map("map.map");
}
