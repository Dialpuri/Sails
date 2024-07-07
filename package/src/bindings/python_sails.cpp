//
// Created by Jordan Dialpuri on 16/02/2024.
//
#include <nanobind/nanobind.h>

#include "../include/sails-json.h"
#include "../include/sails-sequence.h"
#include "../include/sails-glycan.h"
#include "../include/sails-topology.h"

#include "gemmi/model.hpp" // for Structure
#include "gemmi/mmread.hpp" // for read_structure
#include "gemmi/resinfo.hpp" // for find_tabulated_residue

namespace nb = nanobind;
using namespace nb::literals;


void run() {
    Sails::ResidueDatabase database;
    Sails::JSONLoader loader = {"package/data/data.json", database};

    const std::string path = "package/models/5fji/5fji.cif";
    gemmi::Structure structure = gemmi::read_structure_file(path);

    Sails::Glycosites glycosites = Sails::find_n_glycosylation_sites(structure);
    std::string dot_file;

    for (auto& glycosite: glycosites) {
        auto glycan = Sails::find_glycan_topology(structure, glycosite, database);

        if (glycan.has_value()) {
            dot_file += glycan->get_dot_string(structure, database);
        }
    }
    std::ofstream of("package/graphs/all-glycans.dot");
    of << dot_file << std::endl;
    of.close();

    // generate tree structure of glycans
    // check glycosite residues
    // check donor atoms
    // check for attached sugar with neighbour search
    // check donors of attached sugar for other sugars
    // tree structure
    //
}

NB_MODULE(sails_module, m) {
    m.def("run", &run);
}