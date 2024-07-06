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

    for (const auto& glycosite: glycosites) {
        auto glycan = Sails::find_glycan_topology(structure, glycosite);
    }

    std::cout << glycosites.size() << " glycosites found." << std::endl;

    Sails::Sugar sugar1 = {"C1", 101};
    Sails::Sugar sugar2 = {"C2", 102};
    Sails::Sugar sugar3 = {"C3", 103};
    Sails::Sugar sugar4 = {"C4", 104};

    Sails::Glycan glycan;
    glycan.addLinkage(sugar1, sugar2);
    glycan.addLinkage(sugar2, sugar3);
    glycan.addLinkage(sugar1, sugar4);



//    glycan.print_list();

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