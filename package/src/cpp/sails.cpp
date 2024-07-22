//
// Created by Jordan Dialpuri on 06/07/2024.
//

#include "../include/sails-density.h"
#include "../include/sails-json.h"
#include "../include/sails-sequence.h"
#include "../include/sails-glycan.h"
#include "../include/sails-topology.h"
#include "../include/sails-linkage.h"
#include "../include/sails-cif.h"


#include "gemmi/model.hpp" // for Structure
#include "gemmi/mmread.hpp" // for read_structure
#include "gemmi/resinfo.hpp" // for find_tabulated_residue

#include <chrono>
#include <iostream>

void display_progress_bar(size_t total, size_t &progress_count) {
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

void remove_erroneous_sugars(gemmi::Structure *structure, Sails::Density *density, Sails::Glycan *glycan) {
    float rscc_threshold = 0.4;
    std::vector<Sails::Sugar *> to_remove;
    for (auto &sugar: *glycan) {
        auto residue = Sails::Utils::get_residue_from_glycosite(sugar.second->site, structure);
        auto sugar_result = glycan->find_previous_sugar(sugar.second.get());
        if (!sugar_result) continue;
        auto previous_residue =
                Sails::Utils::get_residue_from_glycosite(sugar_result.value()->site, structure);

        if (residue.name == "ASN") { continue; } // don't remove ASN
        if (residue.name == "TRP") { continue; } // don't remove ASN

        if (float rscc = round(10 * density->rscc_score(residue))/10; rscc < rscc_threshold) {
            to_remove.emplace_back(sugar.second.get()); // add pointer to
            std::cout << "Removing " << Sails::Utils::format_residue_key(&residue) << " because of low RSCC =" << rscc
                    << std::endl;
            continue;
        }

        // cases where sugar is clashing into protein density
        if (float diff_score = density->difference_density_score(residue); diff_score > 1.1) {
            std::cout << "Removing " << Sails::Utils::format_residue_key(&previous_residue) << "--" <<
                    Sails::Utils::format_residue_key(&residue) << " because of high DDS = " << diff_score <<
                    std::endl;
            to_remove.emplace_back(sugar.second.get());
        }
    }

    // sort removal in decsending order so removed indices don't cause later array overflow
    std::sort(to_remove.begin(), to_remove.end(), [](const Sails::Sugar *a, const Sails::Sugar *b) {
        return !(a->site < b->site);
    });

    for (auto &sugar: to_remove) {
        glycan->remove_sugar(sugar);
    }
}

std::vector<Sails::LinkRecord> generate_link_records(gemmi::Structure *structure, Sails::Glycosites *glycosites, Sails::Topology *topology) {
    std::vector<Sails::LinkRecord> links;
    int link_id = 0;
    for (auto &glycosite: *glycosites) {
        Sails::Glycan glycan = topology->find_glycan_topology(glycosite);
        std::vector<Sails::Linkage> list = glycan.linkage_list;

        for (auto l : list) {
            std::string id = "covalent" + std::to_string(link_id++);
            gemmi::Chain c1 = Sails::Utils::get_chain_from_glycosite(l.donor_sugar->site, structure);
            gemmi::Chain c2 = Sails::Utils::get_chain_from_glycosite(l.acceptor_sugar->site, structure);

            gemmi::Residue r1 = Sails::Utils::get_residue_from_glycosite(l.donor_sugar->site, structure);
            gemmi::Residue r2 = Sails::Utils::get_residue_from_glycosite(l.acceptor_sugar->site, structure);

            gemmi::Atom* a1 = &r1.sole_atom(l.donor_atom);
            gemmi::Atom* a2 = &r2.sole_atom(l.acceptor_atom);

            Sails::LinkRecord x = {id, c1, c2, r1, r2, *a1, *a2};
            links.push_back(x);
        }
    }
    return links;
}

void run() {

    Sails::JSONLoader loader = {"package/data/data.json"};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    // N-glycosylation
    const std::string path = "package/models/5fji/5fji_deglycosylated.pdb";
    // const std::string path = "package/models/5fji/5fji.cif";
    const std::string mtz_path = "package/models/5fji/5fji.mtz";

    // C-Mannosylation Tests
    // const std::string path = "package/models/6plh/6plh_deglycosylated.cif";
    // const std::string mtz_path = "package/models/6plh/6plh.mtz";

    gemmi::Structure structure = gemmi::read_structure_file(path);
    gemmi::Structure original_structure = structure;

    Sails::Glycosites glycosites = Sails::find_n_glycosylation_sites(structure);
    // Sails::Glycosites glycosites = Sails::find_c_glycosylation_sites(structure);
    // auto site = Sails::find_site(structure, "A", "ASN", 323);

    // auto site = Sails::find_site(structure, "B", "ASN", 524);
    // if (!site.has_value()) throw std::runtime_error("Site not found");
    // Sails::Glycosites glycosites = {site.value()};

    Sails::Topology topology = {&structure, residue_database};

    Sails::Density density = Sails::Density(mtz_path);
    density.load_hkl("FP", "SIGFP");
    density.recalculate_map(structure);
    density.calculate_alt_map(original_structure);

    structure.cell = density.m_mtz.cell;
    structure.spacegroup_hm = density.m_mtz.spacegroup_name;

    Sails::Model model = {&structure, linkage_database, residue_database};

    int cycles = 1;
    size_t progress_count = 0;

    for (int i = 1; i <= cycles; i++) {
        std::cout << "Cycle #" << i << std::endl;
        // display_progress_bar(cycles, progress_count);

        for (auto &glycosite: glycosites) {
            // display the progress bar

            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            if (glycan.empty()) {
                continue;
            }

            // find terminal sugars
            // auto residue = Sails::Utils::get_residue_from_glycosite(glycosite, &structure);
            Sails::Glycan new_glycan = model.extend(glycan, glycosite, density);
            topology.set_structure(model.get_structure());
        }

        // recalculate the map once all sugars have been added
        // then go back through all the glycans and remove the sugars which are bad
        // MAP RECALCULATION
        density.recalculate_map(structure);
        density.calculate_alt_map(original_structure);

        // remove erroneous sugars
        for (auto &glycosite: glycosites) {
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            if (glycan.empty()) { continue; }

            remove_erroneous_sugars(&structure, &density, &glycan);
        }
    }


    // LINKS - NOT FINISHED
    gemmi::cif::Loop loop;
    loop.tags = Sails::LinkRecord::tags();
    std::vector<Sails::LinkRecord> links = generate_link_records(&structure, &glycosites, &topology);
    model.save("structure.cif", links);
    density.write_mtz("reflections.mtz");
}

int main() {
    run();
}

