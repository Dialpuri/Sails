//
// Created by Jordan Dialpuri on 13/08/2024.
//

#include "../../include/density/sails-density.h"
#include "../../include/density/sails-em-density.h"

Sails::EMDensity::EMDensity(gemmi::Grid<> &grid) {
    m_grid = grid;
}
