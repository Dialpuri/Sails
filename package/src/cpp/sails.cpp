//
// Created by Jordan Dialpuri on 06/07/2024.
//

#include "../include/sails-json.h"
#include "../include/sails-sequence.h"
#include "../include/sails-glycan.h"
#include "../include/sails-topology.h"
#include "../include/sails-linkage.h"

#include "gemmi/model.hpp" // for Structure
#include "gemmi/mmread.hpp" // for read_structure
#include "gemmi/resinfo.hpp" // for find_tabulated_residue


void run() {
    Sails::JSONLoader loader = {"package/data/data.json"};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    const std::string path = "package/models/5fji/5fji.cif";
    gemmi::Structure structure = gemmi::read_structure_file(path);

    Sails::Glycosites glycosites = Sails::find_n_glycosylation_sites(structure);
    Sails::Topology topology = {structure, residue_database};

    Sails::Model model = {structure, linkage_database,residue_database};
    for (auto& glycosite: glycosites) {
        auto glycan = topology.find_glycan_topology(glycosite);
        auto r = Sails::Utils::get_residue_from_glycosite(glycosite, structure);
        if (glycan.empty()) {continue;}

        // remove sugars with a lot of red difference density nearby -> indicates its wrong and then store those so we dont try rebuilding there...

        // find terminal sugars
        auto residue = Sails::Utils::get_residue_from_glycosite(glycosite, structure);
        auto new_glycan = model.extend(glycan, residue.seqid.num.value);
    }

}

int main() {
    run();
}