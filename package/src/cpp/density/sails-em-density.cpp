//
// Created by Jordan Dialpuri on 13/08/2024.
//

#include "../../include/density/sails-density.h"
#include "../../include/density/sails-em-density.h"

Sails::EMDensity::EMDensity(gemmi::Grid<> &grid, float resolution) {
    m_grid = grid;
    m_resolution = resolution;
}

gemmi::Grid<> Sails::EMDensity::calculate_density_for_box(gemmi::Residue &residue,
    gemmi::Box<gemmi::Position> &box) const {
    gemmi::DensityCalculator<gemmi::C4322<float>, float> density_calculator;

    gemmi::Position size = box.get_size();
    gemmi::UnitCell dummy_cell = {size.x, size.y, size.z, 90, 90, 90};
    density_calculator.grid.unit_cell = dummy_cell;
    density_calculator.grid.nu = size.x;
    density_calculator.grid.nv = size.y;
    density_calculator.grid.nw = size.z;
    density_calculator.grid.spacegroup = get_work_grid()->spacegroup;
    density_calculator.grid.axis_order = get_work_grid()->axis_order;

    density_calculator.d_min = 1;
    density_calculator.initialize_grid();
    for (auto &atom: residue.atoms) {
        density_calculator.add_atom_density_to_grid(atom);
    }
    density_calculator.grid.symmetrize_sum();
    return density_calculator.grid;
}

gemmi::Grid<> Sails::EMDensity::calculate_density_for_grid(gemmi::Residue &residue) const {
    gemmi::DensityCalculator<gemmi::C4322<float>, float> density_calculator;

    density_calculator.grid.copy_metadata_from(*get_best_grid());
    density_calculator.grid.spacing[0] = get_best_grid()->spacing[0];
    density_calculator.grid.spacing[1] = get_best_grid()->spacing[1];
    density_calculator.grid.spacing[2] = get_best_grid()->spacing[2];

    density_calculator.d_min = get_resolution();
    density_calculator.initialize_grid();
    for (auto &atom: residue.atoms) {
        density_calculator.add_atom_density_to_grid(atom);
    }
    density_calculator.grid.symmetrize_sum();
    auto x =  density_calculator.grid;
    return std::move(x);
}

gemmi::Grid<> Sails::EMDensity::calculate_density_for_structure(gemmi::Structure &structure) const {
     gemmi::DensityCalculator<gemmi::C4322<float>, float> density_calculator;

     density_calculator.grid.copy_metadata_from(*get_best_grid());
     density_calculator.grid.spacing[0] = get_best_grid()->spacing[0];
     density_calculator.grid.spacing[1] = get_best_grid()->spacing[1];
     density_calculator.grid.spacing[2] = get_best_grid()->spacing[2];

     density_calculator.d_min = get_resolution();
     density_calculator.initialize_grid();
     density_calculator.add_model_density_to_grid(structure.models[0]);
     density_calculator.grid.symmetrize_sum();
     auto x =  density_calculator.grid;
     return std::move(x);
}
