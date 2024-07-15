//
// Created by Jordan Dialpuri on 06/07/2024.
//

#include "../include/sails-density.h"
#include "../include/sails-json.h"
#include "../include/sails-sequence.h"
#include "../include/sails-glycan.h"
#include "../include/sails-topology.h"
#include "../include/sails-linkage.h"

#include "gemmi/model.hpp" // for Structure
#include "gemmi/mmread.hpp" // for read_structure
#include "gemmi/resinfo.hpp" // for find_tabulated_residue

#include <chrono>
#include <iostream>

void display_progress_bar(size_t total, size_t& progress_count) {
    std::cout << "[";
    int position = progress_count * 60 / total;
    for (int i = 0; i < 60; ++i) {
        if (i < position) std::cout << "=";
        else if (i == position) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int((progress_count * 100.0) / total) << " %\r";
    std::cout.flush();
    ++progress_count;
}

void run() {
    Sails::JSONLoader loader = {"package/data/data.json"};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    const std::string path = "package/models/5fji/5fji_deglycosylated.pdb";
    // const std::string path = "package/models/5fji/5fji.cif";

    gemmi::Structure structure = gemmi::read_structure_file(path);

    Sails::Glycosites glycosites = Sails::find_n_glycosylation_sites(structure);
    Sails::Topology topology = {structure, residue_database};

    const std::string mtz_path = "package/models/5fji/5fji.mtz";
    Sails::Density density = Sails::Density(mtz_path, "FWT", "PHWT");

    Sails::Model model = {structure, linkage_database, residue_database};

    size_t total = glycosites.size();
    int cycles = 3;

    for (int i = 0; i < cycles; i++) {
        std::cout << "Cycle #" << i << std::endl;
        size_t progress_count = 0;

        for (auto &glycosite: glycosites) {
            // display the progress bar
            // display_progress_bar(total, progress_count);

            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            if (glycan.empty()) { continue; }

            // find terminal sugars
            auto residue = Sails::Utils::get_residue_from_glycosite(glycosite, structure);
            Sails::Glycan new_glycan = model.extend(glycan, residue.seqid.num.value, density);

            // set the topology to have the updated structure
            gemmi::Structure current_structure = model.get_structure();
            topology.set_structure(current_structure);
        }

        structure = topology.get_structure(); // update to the latest structure

        // recalculate the map once all sugars have been added
        // then go back through all the glycans and remove the sugars which are bad

        // MAP RECALCULATION

        // remove erroneous sugars
        for (auto &glycosite: glycosites) {
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            if (glycan.empty()) { continue; }

            float rscc_threshold = 0.3;
            std::vector<Sails::Sugar*> to_remove;
            for (auto& sugar: glycan) {
                auto residue = Sails::Utils::get_residue_from_glycosite(sugar.second->site, structure);
                if (float rscc = density.rscc_score(residue); rscc < rscc_threshold) {
                    to_remove.emplace_back(sugar.second.get()); // add pointer to
                }
            }

            std::vector<Sails::Sugar*> unbuildable_terminals;
            for (auto& sugar: to_remove) {
                glycan.remove_sugar(sugar);
            }

            gemmi::Structure current_structure = glycan.get_structure();
            topology.set_structure(current_structure); // set the topology to have the updated structure
        }
    }

    model.save_pdb("structure.pdb");
}

int main() {
    run();
}
