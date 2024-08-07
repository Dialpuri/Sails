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
#include "../include/sails-telemetry.h"
#include "../include/snfg/sails-snfg.h"

#include "gemmi/model.hpp" // for Structure
#include "gemmi/mmread.hpp" // for read_structure
#include "gemmi/resinfo.hpp" // for find_tabulated_residue

#include <chrono>
#include <iostream>
#include <src/include/sails-gemmi-bindings.h>


void print_rejection_dds(const Sails::Glycosite& s1, const Sails::Glycosite& s2, gemmi::Structure* structure, float score) {
    std::cout << "Removing " << Sails::Utils::format_residue_from_site(s1, structure) << "--"
    << Sails::Utils::format_residue_from_site(s2, structure) << " because of high DDS = " << score <<std::endl;
}

void print_removal_rscc(const gemmi::Residue& residue, float rscc) {
    std::cout << "Removing " << Sails::Utils::format_residue_key(&residue) << " because of low RSCC =" << rscc << std::endl;
}

void remove_erroneous_sugars(gemmi::Structure *structure, Sails::Density *density, Sails::Glycan *glycan, bool debug) {
    constexpr float rscc_threshold = 0.5;
    constexpr float dds_threshold = 1.1;

    std::vector<Sails::Sugar *> to_remove;
    for (const auto &[fst, snd]: *glycan) {
        gemmi::Residue residue = Sails::Utils::get_residue_from_glycosite(snd->site, structure);

        std::optional<Sails::Sugar *> sugar_result = glycan->find_previous_sugar(snd.get());
        if (!sugar_result.has_value()) continue; // if there is nothing previous, it must be a protein residue

        gemmi::Residue previous_residue = Sails::Utils::get_residue_from_glycosite(
            sugar_result.value()->site, structure);

        // if (residue.name == "ASN") { continue; } // don't remove ASN
        // if (residue.name == "TRP") { continue; } // don't remove TRP

        // remove cases with low rscc
        if (const float rscc = density->rscc_score(residue); rscc < rscc_threshold) {
            to_remove.emplace_back(snd.get()); // add pointer to
            if (debug) print_removal_rscc(residue, rscc);
            continue;
        }

        // remove cases with high difference density score
        if (const float diff_score = density->difference_density_score(residue); diff_score > dds_threshold) {
            if (debug) print_rejection_dds(sugar_result.value()->site, fst, structure, diff_score);
            to_remove.emplace_back(snd.get());
        }
    }

    // add linked sugars to removal list
    for (auto &sugar: to_remove) {
        std::vector<Sails::Sugar *> additional_sugars;
        for (auto &linked_sugar: glycan->adjacency_list[sugar]) {
            // check that the linked sugar is not already in the removal list
            if (std::find(to_remove.begin(), to_remove.end(), linked_sugar) != to_remove.end()) continue;

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

Sails::Output run_cycle(Sails::Glycosites& glycosites, gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, std::string &
                        resource_dir, bool verbose) {

    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    gemmi::Structure original_structure = structure;
    gemmi::Mtz mtz = form_gemmi_mtz(sails_mtz);
    check_spacegroup(&mtz, &structure); // check to ensure the MTZ has a spacegroup

    Sails::Topology topology = {&structure, residue_database};
    Sails::SNFG snfg = Sails::SNFG(&structure, &residue_database);

    Sails::Density density = Sails::Density(mtz);
    density.load_hkl("FP", "SIGFP");
    density.recalculate_map(structure);
    density.calculate_po_pc_map(original_structure);

    structure.cell = density.m_mtz.cell;
    structure.spacegroup_hm = density.m_mtz.spacegroup_name;

    Sails::Model model = {&structure, linkage_database, residue_database};
    model.set_special_monomer_dir(resource_dir);

    Sails::Telemetry telemetry = Sails::Telemetry("");

    for (int i = 1; i <= cycles; i++) {
        if (!verbose) std::cout << "\rCycle #" << i;
        std::cout << std::flush;
        if (verbose) std::cout << "\rCycle #" << i << std::endl;

        for (auto &glycosite: glycosites) {
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            // if (glycan.empty()) { continue; }

            // find terminal sugars
            Sails::Glycan new_glycan = model.extend(glycan, glycosite, density, verbose);

            std::set<Sails::Glycosite> differences = new_glycan - glycan;
            telemetry << differences;

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
            Sails::Glycan old_glycan = glycan;
            remove_erroneous_sugars(&structure, &density, &glycan, verbose);

            topology.set_structure(&structure); // need to update neighbor search after removing n residues
            Sails::Glycan new_glycan = topology.find_glycan_topology(glycosite);

            std::set<Sails::Glycosite> differences = old_glycan - new_glycan;
            telemetry >> differences;

            std::string snfg_string = snfg.create_snfg(new_glycan, glycosite);
            std::string glycosite_key = Sails::Utils::format_residue_from_site(glycosite, &structure);
            telemetry.save_snfg(i, glycosite_key, snfg_string);
        }

        telemetry.save_state(i);
    }

    std::cout << std::endl;

    // add links and write files
    std::vector<Sails::LinkRecord> links = generate_link_records(&structure, &glycosites, &topology);

    Sails::MTZ output_mtz = Sails::form_sails_mtz(density.m_mtz, "FP", "SIGFP");
    std::string log_string = telemetry.format_log(&structure, &density, false).value();

    Sails::Telemetry::SNFGCycleData snfgs = telemetry.get_snfgs();
    return {
        *model.get_structure(),
        output_mtz,
        log_string,
        snfgs
    };
}

Sails::Output n_glycosylate(gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, std::string &resource_dir,
                            bool verbose) {
    auto glycosites = Sails::find_n_glycosylation_sites(structure);
    return run_cycle(glycosites, structure, sails_mtz, cycles, resource_dir, verbose);
}

Sails::Output c_glycosylate(gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, std::string &resource_dir,
                            bool verbose) {
    auto glycosites = Sails::find_c_glycosylation_sites(structure);
    return run_cycle(glycosites, structure, sails_mtz, cycles, resource_dir, verbose);
}

std::string get_snfg(std::string chain, int seqid, gemmi::Structure& structure, std::string& resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();

    Sails::Topology topology = {&structure, residue_database};
    Sails::SNFG snfg = Sails::SNFG(&structure, &residue_database);

    std::optional<Sails::Glycosite> potential_glycosite = Sails::find_site(structure, chain, seqid);
    if (!potential_glycosite.has_value()) throw std::runtime_error("Could not find specified site");
    Sails::Glycosite glycosite = potential_glycosite.value();

    Sails::Glycan glycan = topology.find_glycan_topology(glycosite);

    return snfg.create_snfg(glycan, glycosite);
}

std::map<std::string, std::string> get_all_snfgs(gemmi::Structure& structure, std::string& resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();

    Sails::Topology topology = {&structure, residue_database};
    Sails::SNFG snfg = Sails::SNFG(&structure, &residue_database);

    std::map<std::string, std::string> snfg_map;
    Sails::Glycosites n_glycosites = Sails::find_n_glycosylation_sites(structure);
    for (auto& site: n_glycosites) {
        Sails::Glycan glycan = topology.find_glycan_topology(site);
        if (glycan.empty()) continue;
        std::string key = Sails::Utils::format_residue_from_site(site, &structure);
        snfg_map[key] = snfg.create_snfg(glycan, site);
    }

    Sails::Glycosites c_glycosites = Sails::find_c_glycosylation_sites(structure);
    for (auto& site: c_glycosites) {
        Sails::Glycan glycan = topology.find_glycan_topology(site);
        if (glycan.empty()) continue;
        std::string key = Sails::Utils::format_residue_from_site(site, &structure);
        snfg_map[key] = snfg.create_snfg(glycan, site);
    }

    return snfg_map;
}

void test() {
    const std::string path = "testing/test_data/5fji/5FJI.cif";
    const std::string mtz_path = "testing/test_data/5fji/5fji.mtz";

    gemmi::Structure structure = gemmi::read_structure_file(path);

    std::string data_file = "package/src/sails/data/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();

    auto snfg = Sails::SNFG(&structure, &residue_database);
    Sails::Topology topology = {&structure, residue_database};
    auto glycosites = Sails::find_n_glycosylation_sites(structure);

    for (auto& site: glycosites) {
        // auto site = glycosites[4];
        auto glycan = topology.find_glycan_topology(site);
        if (glycan.empty()) continue;
        std::string snfg_path = "snfgs/" + Sails::Utils::format_residue_from_site(site, &structure) + ".svg";
        std::ofstream f(snfg_path);
        f << snfg.create_snfg(glycan, site);
    f.close();
    }
}

// testbed
int main() {
    const std::string path = "testing/test_data/5fji/5fji.cif";
    const std::string mtz_path = "testing/test_data/5fji/5fji.mtz";

    gemmi::Structure structure = gemmi::read_structure_file(path);
    gemmi::Mtz mtz = gemmi::read_mtz_file(mtz_path);
    Sails::MTZ sails_mtz = Sails::form_sails_mtz(mtz, "FP", "SIGFP");
    std::string data_file = "package/src/sails/data/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();

    auto snfg = Sails::SNFG(&structure, &residue_database);
    Sails::Topology topology = {&structure, residue_database};
    auto glycosites = Sails::find_n_glycosylation_sites(structure);

    for (auto &site: glycosites) {
        auto glycan = topology.find_glycan_topology(site);
        snfg.create_snfg(glycan, site);
    }
}
