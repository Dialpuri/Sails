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
#include <src/include/sails-gemmi-bindings.h>

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

void remove_erroneous_sugars(gemmi::Structure *structure, Sails::Density *density, Sails::Glycan *glycan, bool debug) {
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

        if (float rscc = round(10 * density->rscc_score(residue)) / 10; rscc < rscc_threshold) {
            to_remove.emplace_back(sugar.second.get()); // add pointer to
            if (debug)
                std::cout << "Removing " << Sails::Utils::format_residue_key(&residue) << " because of low RSCC =" <<
                        rscc << std::endl;
            continue;
        }

        // cases where sugar is clashing into protein density
        if (float diff_score = density->difference_density_score(residue); diff_score > 1.1) {
            if (debug)
                std::cout << "Removing " << Sails::Utils::format_residue_from_site(
                            sugar_result.value()->site, structure) << "--" <<
                        Sails::Utils::format_residue_from_site(sugar.first, structure) << " because of high DDS = " <<
                        diff_score <<
                        std::endl;
            to_remove.emplace_back(sugar.second.get());
        }
    }

    // add linked sugars to removal list
    for (auto &sugar: to_remove) {
        std::vector<Sails::Sugar *> additional_sugars;
        for (auto &linked_sugar: glycan->adjacency_list[sugar]) {
            additional_sugars.emplace_back(linked_sugar);
        }
        to_remove.insert(to_remove.end(), additional_sugars.begin(), additional_sugars.end());
    }

    // sort removal in decsending order so removed indices don't cause later array overflow
    std::sort(to_remove.begin(), to_remove.end(), [](const Sails::Sugar *a, const Sails::Sugar *b) {
        return !(a->site < b->site);
    });

    for (auto &sugar: to_remove) {
        glycan->remove_sugar(sugar);
    }
}

Sails::Glycan get_glycan_topology(gemmi::Structure &structure, Sails::Glycosite &glycosite) {
    Sails::JSONLoader loader = {"package/data/data.json"};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::Topology topology = {&structure, residue_database};
    return topology.find_glycan_topology(glycosite);
}

void check_spacegroup(gemmi::Mtz* mtz, gemmi::Structure* structure) {
    if (!mtz->spacegroup_name.empty() && !structure->spacegroup_hm.empty()) return;
    if (mtz->spacegroup_name.empty() && structure->spacegroup_hm.empty()) throw std::runtime_error("No spacegroup information in MTZ or Structure");
    if (mtz->spacegroup_name.empty()) mtz->spacegroup_name = structure->spacegroup_hm;
}

Sails::Output run_cycle(Sails::Glycosites& glycosites, gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, bool verbose) {
    Sails::JSONLoader loader = {"/Users/dialpuri/Development/sails/package/data/data.json"};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    gemmi::Structure original_structure = structure;
    gemmi::Mtz mtz = form_gemmi_mtz(sails_mtz);
    check_spacegroup(&mtz, &structure); // check to ensure the MTZ has a spacegroup

    Sails::Topology topology = {&structure, residue_database};

    Sails::Density density = Sails::Density(mtz);
    density.load_hkl("FP", "SIGFP");
    density.recalculate_map(structure);
    density.calculate_po_pc_map(original_structure);

    structure.cell = density.m_mtz.cell;
    structure.spacegroup_hm = density.m_mtz.spacegroup_name;

    Sails::Model model = {&structure, linkage_database, residue_database};

    size_t progress_count = 0;

    for (int i = 1; i <= cycles; i++) {
        if (!verbose) std::cout << "\rCycle #" << i;
        std::cout << std::flush;
        if (verbose) std::cout << "\rCycle #" << i << std::endl;
        // display_progress_bar(cycles, progress_count);

        for (auto &glycosite: glycosites) {
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            if (glycan.empty()) { continue; }

            // find terminal sugars
            Sails::Glycan new_glycan = model.extend(glycan, glycosite, density, verbose);
            topology.set_structure(model.get_structure());
        }

        // recalculate maps
        density.recalculate_map(structure);
        density.calculate_po_pc_map(original_structure);

        // remove erroneous sugars
        for (auto &glycosite: glycosites) {
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            if (glycan.empty()) { continue; }

            // std::cout << "Attempting removal at " << Sails::Utils::format_residue_from_site(glycosite, &structure) << std::endl;
            remove_erroneous_sugars(&structure, &density, &glycan, verbose);
        }
    }

    std::cout << std::endl;
    // add links and write files
    std::vector<Sails::LinkRecord> links = generate_link_records(&structure, &glycosites, &topology);

    Sails::MTZ output_mtz = Sails::form_sails_mtz(density.m_mtz, "FP", "SIGFP");
    return {
        *model.get_structure(),
        output_mtz
    };
}

Sails::Output n_glycosylate(gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, bool verbose) {
    auto glycosites = Sails::find_n_glycosylation_sites(structure);
    return run_cycle(glycosites, structure, sails_mtz, cycles, verbose);
}

Sails::Output c_glycosylate(gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, bool verbose) {
    auto glycosites = Sails::find_c_glycosylation_sites(structure);
    return run_cycle(glycosites, structure, sails_mtz, cycles, verbose);
}

int main() {
    // const std::string path = "/Users/dialpuri/Development/sails/testing/test_data/8dvl/8DVL_deglycosylated.cif";
    // const std::string mtz_path = "/Users/dialpuri/Development/sails/testing/test_data/8dvl/8DVL.mtz";

    const std::string path = "/Users/dialpuri/Development/sails/testing/test_data/5fji/5FJI_deglycosylated.cif";
    const std::string mtz_path = "/Users/dialpuri/Development/sails/testing/test_data/5fji/5FJI.mtz";
    //
    // const std::string path = "/Users/dialpuri/Development/sails/testing/test_data/3sku/3SKU_deglycosylated.cif";
    // const std::string mtz_path = "/Users/dialpuri/Development/sails/testing/test_data/3sku/3SKU.mtz";
    int cycles = 3;

    gemmi::Structure structure = gemmi::read_structure_file(path);
    gemmi::Mtz mtz = gemmi::read_mtz_file(mtz_path);
    Sails::MTZ sails_mtz = Sails::form_sails_mtz(mtz, "FP", "SIGFP");
    auto output = n_glycosylate(structure, sails_mtz, cycles, true);

    std::string output_path = "structure.cif";
    std::string output_mtz_path = "reflections.mtz";

    std::ofstream os(output_path);
    gemmi::cif::Document document = make_mmcif_document(output.structure);
    write_cif_to_stream(os, document);
    os.close();

    form_gemmi_mtz(output.mtz).write_to_file(output_mtz_path);
}
