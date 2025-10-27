//
// Created by Jordan Dialpuri on 06/07/2024.
//

#include "../include/density/sails-density.h"
#include "../include/density/sails-xtal-density.h"
#include "../include/density/sails-em-density.h"
#include "../include/sails-json.h"
#include "../include/sails-sequence.h"
#include "../include/sails-glycan.h"
#include "../include/sails-topology.h"
#include "../include/sails-linkage.h"
#include "../include/sails-cif.h"
#include "../include/sails-telemetry.h"
#include "../include/sails-wurcs.h"
#include "../include/snfg/sails-snfg.h"
#include <src/include/sails-gemmi-bindings.h>
#include <src/include/sails-solvent.h>

#include "gemmi/model.hpp" // for Structure
#include "gemmi/mmread.hpp" // for read_structure
#include "gemmi/resinfo.hpp" // for find_tabulated_residue
#include "gemmi/ccp4.hpp" // for find_tabulated_residue

#include <chrono>
#include <iostream>
#include <src/include/sails-morph.h>

#include "src/include/sails-predictions.h"
#include "src/include/sails-score.h"


void print_rejection_dds(const Sails::Glycosite& s1, const Sails::Glycosite& s2, gemmi::Structure* structure) {
    std::cout << "Removing " << Sails::Utils::format_residue_from_site(s1, structure) << "--"
    << Sails::Utils::format_residue_from_site(s2, structure) << " because of negative difference density " << std::endl;
}

void print_removal_rscc(const Sails::Glycosite &site, float rscc, gemmi::Structure *structure) {
    std::cout << "Removing " << Sails::Utils::format_residue_from_site(site, structure) << " because of low RSCC =" << rscc << std::endl;
}

void print_rscc(const Sails::Glycosite &site, float rscc, gemmi::Structure *structure) {
    std::cout << Sails::Utils::format_residue_from_site(site, structure) << " - RSCC = " << rscc << std::endl;
}

void print_removal_clash(const Sails::Glycosite &site, float rscc, gemmi::Structure *structure) {
    std::cout << "Removing " << Sails::Utils::format_residue_from_site(site, structure) << " because of clashes (Clash score = " << rscc << ")" << std::endl;
}

void print_dds(const Sails::Glycosite &site, float dds, gemmi::Structure *structure) {
    std::cout << Sails::Utils::format_residue_from_site(site, structure) << " - DDS = " << dds << std::endl;
}

void remove_erroneous_sugars(gemmi::Structure *structure, Sails::Density *density, Sails::Glycan *glycan, bool strict,
                             bool debug, Sails::ResidueDatabase &residue_database) {
    const float rscc_threshold = strict ? 0.65: 0.5;
    const float dds_threshold = strict ? 1.0: 1.1;

    // const std::pair<float, float> difference_density_stats = density->calculate_map_statistics(density->get_difference_grid());
    std::map<Sails::Glycosite, double> rsccs = Sails::Score::calculate_rsccs(density, structure, residue_database);
    std::map<Sails::Glycosite, double> qscores = Sails::Score::calculate_qscores(density, structure, residue_database);

    std::vector<Sails::Sugar *> to_remove;
    for (const auto &[fst, snd]: *glycan) {
        gemmi::Residue residue = Sails::Utils::get_residue_from_glycosite(snd->site, structure);

        std::optional<Sails::Sugar *> sugar_result = glycan->find_previous_sugar(snd.get());
        if (!sugar_result.has_value()) continue; // if there is nothing previous, it must be a protein residue

        if (residue.name == "FUC") {
            double clash_score = Sails::Score::calculate_clash_score(&residue, structure);
            if (clash_score > 2) {
                print_removal_clash(snd->site, clash_score, structure) ;
                to_remove.push_back(snd.get());
                continue;
            }
        }

        gemmi::Residue previous_residue = Sails::Utils::get_residue_from_glycosite(
            sugar_result.value()->site, structure);

        snd->site.atom_idx = 0; // set atom index to 0 so can be used in comparisons on the residue level

        // remove cases with low rscc
        if (rsccs.count(snd->site) != 0) {
            const double rscc = rsccs.at(snd->site);
            print_rscc(snd->site, rscc, structure);
            if (rscc < rscc_threshold) {
                to_remove.emplace_back(snd.get()); // add pointer to remove
                if (debug) print_removal_rscc(snd->site, rscc, structure);
            }
        } else {
            std::cout << Sails::Utils::format_site_key(fst) << " | " << Sails::Utils::format_site_key(snd->site) << std::endl;
            throw std::runtime_error("Glycosite was not found in the RSCC calculation" + Sails::Utils::format_residue_from_site(snd->site, structure));
        }

        // remove cases with high difference density score
        // const int no_atoms_in_negative_density = density->check_difference_density(residue, difference_density_stats);
        // // std::cout << Sails::Utils::format_residue_from_site(fst, structure) << " " << no_atoms_in_negative_density << std::endl;
        // if (no_atoms_in_negative_density > 4) {
        //     if (debug) print_rejection_dds(sugar_result.value()->site, fst, structure);
        //     to_remove.emplace_back(snd.get());
        // }
        // print_dds(snd->site, diff_score, structure);
        // if ( diff_score > dds_threshold) {

        // }
    }

    // add linked sugars to removal list
    std::set<Sails::Sugar *> additional_sugars;
    for (auto &sugar: to_remove) {
        std::vector<Sails::Sugar*> downstream_sugars = glycan->get_downstream_sugars(sugar);

        for (auto& downstream_sugar: downstream_sugars) {
            if (std::find(to_remove.begin(), to_remove.end(), downstream_sugar) != to_remove.end()) continue;
            additional_sugars.insert(downstream_sugar);
        }
    }
    to_remove.insert(to_remove.end(), additional_sugars.begin(), additional_sugars.end());

    // sort removal in decsending order so removed indices don't cause later array overflow
    std::sort(to_remove.begin(), to_remove.end(), [](const Sails::Sugar *a, const Sails::Sugar *b) {
        return !(a->site < b->site);
    });

    for (const auto &sugar: to_remove) {
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

Sails::Output run_cycle(Sails::Glycosites &glycosites, gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles,
                        std::string &resource_dir, bool strict, bool verbose) {

    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    gemmi::Structure original_structure = structure;
    gemmi::Mtz mtz = form_gemmi_mtz(sails_mtz);
    check_spacegroup(&mtz, &structure); // check to ensure the MTZ has a spacegroup

    Sails::Topology topology = {&structure, residue_database};
    Sails::SNFG snfg = Sails::SNFG(&structure, &residue_database);

    auto density = Sails::XtalDensity(mtz);
    density.recalculate_map(structure);
    density.calculate_po_pc_map(original_structure);

    //
    // gemmi::Grid<> x = *density.get_work_grid();
    // gemmi::Ccp4<> m;
    // m.grid = x ;
    // m.update_ccp4_header();
    // m.write_ccp4_map("wrk.map");

    structure.cell = density.get_mtz()->cell;
    structure.spacegroup_hm = density.get_mtz()->spacegroup_name;

    Sails::Model model = {&structure, linkage_database, residue_database};
    model.set_special_monomer_dir(resource_dir);

    Sails::Telemetry telemetry = Sails::Telemetry("");

    Sails::Glycosites original_glycosites = glycosites;

    for (int i = 1; i <= cycles; i++) {
        if (!verbose) std::cout << "\rCycle #" << i;
        std::cout << std::flush;
        if (verbose) std::cout << "\rCycle #" << i << std::endl;

        if (glycosites.empty()) break;
        std::set<Sails::Glycosite> unmodellable_sites = {};

        for (auto &glycosite: glycosites) {
            // auto c = Sails::Utils::get_chain_from_glycosite(glycosite, &structure);
            // auto r = Sails::Utils::get_residue_from_glycosite(glycosite, &structure);
            // if (c.name != "D" || r.seqid.num.value != 483) continue;
            //
            // std::cout << "Checking " << Sails::Utils::format_residue_from_site(glycosite, &structure) << std::endl;
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            // if (glycan.empty()) { continue; }

            // find terminal sugars
            Sails::Glycan new_glycan = model.extend(glycan, glycosite, density, verbose);

            // if nothing was added, add site to unmodellable list
            if (new_glycan.size() == glycan.size()) {
                std::cout << "Nothing new modelled at site:" << Sails::Utils::format_residue_from_site(glycosite, &structure) << std::endl;
                unmodellable_sites.insert(glycosite);
            }

            std::set<Sails::Glycosite> differences = new_glycan - glycan;
            telemetry << differences;

            topology.set_structure(model.get_structure());
        }

        // recalculate maps
        density.recalculate_map(structure);
        density.calculate_po_pc_map(original_structure);

        // const auto x = density.get_mtz();
        // std::string y = "wrk" + std::to_string(i) + ".mtz";
        // x->write_to_file(y);
        // std::string z = "wrk" + std::to_string(i) + ".cif";
        //
        // Sails::Utils::save_structure_to_file(structure, z);

        // remove erroneous sugars
        for (auto &glycosite: glycosites) {
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
            if (glycan.empty()) { continue; }

            // std::cout << "Attempting removal at " << Sails::Utils::format_residue_from_site(glycosite, &structure) << std::endl;
            Sails::Glycan old_glycan = glycan;
            remove_erroneous_sugars(&structure, &density, &glycan, strict, verbose, residue_database);

            topology.set_structure(&structure); // need to update neighbor search after removing n residues
            Sails::Glycan new_glycan = topology.find_glycan_topology(glycosite);

            if (new_glycan.empty()) {
                unmodellable_sites.insert(glycosite);
                continue;
            }

            new_glycan.renumber();

            std::set<Sails::Glycosite> differences = old_glycan - new_glycan;
            telemetry >> differences;

            std::string snfg_string = snfg.create_snfg(new_glycan, glycosite);
            std::string glycosite_key = Sails::Utils::format_residue_from_site(glycosite, &structure);
            telemetry.save_snfg(i, glycosite_key, snfg_string);
        }

        if (verbose && !unmodellable_sites.empty()) {
            std::cout << "Stopping trials at " << unmodellable_sites.size() << " sites." << std::endl;
            for (const auto& site: unmodellable_sites) {
                std::cout << "\tSite:" << Sails::Utils::format_residue_from_site(site, &structure) << std::endl;
            }
        }

        glycosites.erase(
            std::remove_if(glycosites.begin(), glycosites.end(),[&](const Sails::Glycosite &site) {
                return unmodellable_sites.count(site) > 0;
            }),glycosites.end()
        );

        telemetry.save_state(i);
    }

    std::cout << std::endl;

    model.standardise_residue_names();

    // find and remove any free sugars (likely due to something going wrong)
    std::set<Sails::Glycosite> all_sites = {};
    for (auto &glycosite: original_glycosites) {
        Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
        auto sites = glycan.get_sites();
        all_sites.insert(sites.begin(), sites.end());
    }


    model.remove_free_sites(all_sites);
    topology.set_structure(model.get_structure());

    // add links and write files
    std::vector<Sails::LinkRecord> links = generate_link_records(&structure, &glycosites, &topology);
    Sails::add_links_to_structure(model.get_structure(), links);
    Sails::MTZ output_mtz = Sails::form_sails_mtz(*density.get_mtz(), "FP", "SIGFP");
    std::string log_string = telemetry.format_log(&structure, &density, false).value();

    Sails::Telemetry::SNFGCycleData snfgs = telemetry.get_snfgs();
    return {
        *model.get_structure(),
        output_mtz,
        log_string,
        snfgs
    };
}

Sails::Output run_em_cycle(Sails::Glycosites &glycosites, gemmi::Structure &structure, gemmi::Grid<>& grid, float resolution,  int cycles,
                        std::string &resource_dir, bool strict, bool verbose) {


    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    gemmi::Structure original_structure = structure;

    Sails::Topology topology = {&structure, residue_database};
    Sails::SNFG snfg = Sails::SNFG(&structure, &residue_database);

    auto density = Sails::EMDensity(grid, resolution);

    structure.cell = density.get_mtz()->cell;
    structure.spacegroup_hm = density.get_mtz()->spacegroup_name;

    Sails::Model model = {&structure, linkage_database, residue_database};
    model.set_special_monomer_dir(resource_dir);

    Sails::Telemetry telemetry = Sails::Telemetry("");

    for (int i = 1; i <= cycles; i++) {
        if (!verbose) std::cout << "\rCycle #" << i;
        std::cout << std::flush;
        if (verbose) std::cout << "\rCycle #" << i << std::endl;

        std::cout << "Attempting to model at " << glycosites.size() << " sites." << std::endl;
        for (auto &glycosite: glycosites) {
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);

            // find terminal sugars
            Sails::Glycan new_glycan = model.extend(glycan, glycosite, density, verbose);

            std::set<Sails::Glycosite> differences = new_glycan - glycan;
            telemetry << differences;

            topology.set_structure(&structure);
        }

        // remove erroneous sugars
        std::set<Sails::Glycosite> unmodellable_sites = {};
        for (auto &glycosite: glycosites) {
            Sails::Glycan glycan = topology.find_glycan_topology(glycosite);

            // std::cout << "Attempting removal at " << Sails::Utils::format_residue_from_site(glycosite, &structure) << std::endl;
            Sails::Glycan old_glycan = glycan;
            remove_erroneous_sugars(&structure, &density, &glycan, strict, verbose, residue_database);

            topology.set_structure(&structure); // need to update neighbor search after removing n residues

            Sails::Glycan new_glycan = topology.find_glycan_topology(glycosite);

            if (new_glycan.empty()) {
                unmodellable_sites.insert(glycosite);
                continue;
            }

            new_glycan.renumber();

            std::set<Sails::Glycosite> differences = old_glycan - new_glycan;
            telemetry >> differences;

            std::string snfg_string = snfg.create_snfg(new_glycan, glycosite);
            std::string glycosite_key = Sails::Utils::format_residue_from_site(glycosite, &structure);
            telemetry.save_snfg(i, glycosite_key, snfg_string);
        }

        // sort removal in decsending order so removed indices don't cause later array overflow
        if (verbose && !unmodellable_sites.empty()) {
            std::cout << "Stopping trials at " << unmodellable_sites.size() << " sites." << std::endl;
            for (const auto& site: unmodellable_sites) {
                std::cout << "\tSITE:" << Sails::Utils::format_residue_from_site(site, &structure) << std::endl;
            }
        }

        glycosites.erase(
            std::remove_if(glycosites.begin(), glycosites.end(),[&](const Sails::Glycosite &site) {
                return unmodellable_sites.count(site) > 0;
            }),glycosites.end()
        );

        telemetry.save_state(i);
    }

    std::cout << std::endl;
    model.standardise_residue_names();

    // find and remove any free sugars (likely due to something going wrong)
    std::set<Sails::Glycosite> all_sites = {};
    for (auto &glycosite: glycosites) {
        Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
        auto sites = glycan.get_sites();
        all_sites.insert(sites.begin(), sites.end());
    }

    model.remove_free_sites(all_sites);
    topology.set_structure(model.get_structure());

    // add links and write files
    std::vector<Sails::LinkRecord> links = generate_link_records(&structure, &glycosites, &topology);
    Sails::add_links_to_structure(model.get_structure(), links);
    std::string log_string = telemetry.format_log(&structure, &density, false).value();

    Sails::Telemetry::SNFGCycleData snfgs = telemetry.get_snfgs();
    return {
            *model.get_structure(),
            log_string,
            snfgs
    };
}

Sails::Glycosites identify_predicted_sites(gemmi::Structure &structure, gemmi::Grid<>& glycan_grid, std::string &resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();
    auto predictions = Sails::Predictions(&glycan_grid, linkage_database, residue_database);

    Sails::Glycosites potential_sites = predictions.find_potential_sites(structure, true);
    return potential_sites;
}

Sails::Glycosites identify_predicted_sites(gemmi::Structure &structure, gemmi::Grid<>& glycan_grid, gemmi::Grid<>& protein_grid, bool use_glycan, std::string &resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();
    auto predictions = Sails::Predictions(&glycan_grid, &protein_grid, linkage_database, residue_database);

    Sails::Glycosites potential_sites = predictions.find_potential_sites(structure, use_glycan);
    return potential_sites;
}


// XRAY FUNCTIONS

Sails::Output n_glycosylate(gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, std::string &resource_dir,
                            bool verbose) {
    auto glycosites = Sails::find_n_glycosylation_sites(structure);
    return run_cycle(glycosites, structure, sails_mtz, cycles, resource_dir, true, verbose);
}

Sails::Output c_glycosylate(gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, std::string &resource_dir,
                            bool verbose) {
    auto glycosites = Sails::find_c_glycosylation_sites(structure);
    return run_cycle(glycosites, structure, sails_mtz, cycles, resource_dir, false, verbose);
}

Sails::Output o_mannosylate(gemmi::Structure &structure, Sails::MTZ &sails_mtz, int cycles, std::string &resource_dir,
                            bool verbose) {
    Sails::SolventAccessibility sa = Sails::SolventAccessibility(&structure);
    Sails::SolventAccessibility::SolventAccessibilityMap sa_map = sa.calculate_solvent_accessibility();
    auto glycosites = Sails::find_o_mannosylation_sites(structure, sa_map);
    return run_cycle(glycosites, structure, sails_mtz, cycles, resource_dir, true, verbose);
}

Sails::Output auto_glycosylate(gemmi::Structure &structure, Sails::MTZ &sails_mtz, gemmi::Grid<>& glycan_grid, gemmi::Grid<>& protein_grid, int cycles, std::string &resource_dir,
                            bool verbose) {
    Sails::Glycosites predicted_glycosites = identify_predicted_sites(structure, glycan_grid, protein_grid, true, resource_dir);
    std::cout << "Found " << predicted_glycosites.size() << " potential sites using deep learning models" << std::endl;
    Sails::Glycosites n_glycosites = Sails::find_n_glycosylation_sites(structure);
    Sails::Glycosites c_glycosites = Sails::find_c_glycosylation_sites(structure);

    std::set<Sails::Glycosite> glycosites_set = {predicted_glycosites.begin(), predicted_glycosites.end()};
    glycosites_set.insert(n_glycosites.begin(), n_glycosites.end());
    glycosites_set.insert(c_glycosites.begin(), c_glycosites.end());
    Sails::Glycosites glycosites = {glycosites_set.begin(), glycosites_set.end()};
    int diff = static_cast<int>(glycosites.size()) - static_cast<int>(predicted_glycosites.size());
    std::cout << "Supplemented with " << diff << " sites from the sequence" << std::endl;

    return run_cycle(glycosites, structure, sails_mtz, cycles, resource_dir, false, verbose);
}

Sails::Output glycosylate_site(gemmi::Structure &structure, Sails::MTZ &sails_mtz, std::string& chain, int seqid, int cycles, std::string &resource_dir,
                            bool verbose) {
    std::optional<Sails::Glycosite> potential_site = Sails::find_site(structure, chain, seqid);
    if (!potential_site.has_value()) {
        throw std::runtime_error("Site could not be found");
    }
    Sails::Glycosites glycosites = {potential_site.value()};
    return run_cycle(glycosites, structure, sails_mtz, cycles, resource_dir, false, verbose);
}


// EM FUNCTIONS

Sails::Output n_glycosylate(gemmi::Structure &structure, gemmi::Grid<>& grid, float resolution,  int cycles, std::string &resource_dir,
                            bool verbose) {
    auto glycosites = Sails::find_n_glycosylation_sites(structure);
    return run_em_cycle(glycosites, structure, grid, resolution, cycles, resource_dir, false, verbose);
}

Sails::Output c_glycosylate(gemmi::Structure &structure, gemmi::Grid<>& grid, float resolution,  int cycles, std::string &resource_dir,
                            bool verbose) {
    auto glycosites = Sails::find_c_glycosylation_sites(structure);
    return run_em_cycle(glycosites, structure, grid, resolution, cycles, resource_dir, false, verbose);
}

Sails::Output o_mannosylate(gemmi::Structure &structure, gemmi::Grid<>& grid, float resolution,  int cycles, std::string &resource_dir,
                            bool verbose) {
    Sails::SolventAccessibility sa = Sails::SolventAccessibility(&structure);
    Sails::SolventAccessibility::SolventAccessibilityMap sa_map = sa.calculate_solvent_accessibility();
    auto glycosites = Sails::find_o_mannosylation_sites(structure, sa_map);
    return run_em_cycle(glycosites, structure, grid, resolution, cycles, resource_dir, true, verbose);
}

Sails::Output auto_glycosylate(gemmi::Structure &structure, gemmi::Grid<>& grid, float resolution, gemmi::Grid<>& glycan_grid, gemmi::Grid<>& protein_grid, int cycles, std::string &resource_dir,
                            bool verbose) {
    Sails::Glycosites glycosites = identify_predicted_sites(structure, glycan_grid, protein_grid, true, resource_dir);
    return run_em_cycle(glycosites, structure, grid, resolution, cycles, resource_dir, false, verbose);
}

Sails::Output glycosylate_site(gemmi::Structure &structure, gemmi::Grid<>& grid, float resolution, std::string& chain, int seqid, int cycles, std::string &resource_dir,
                            bool verbose) {
    std::optional<Sails::Glycosite> potential_site = Sails::find_site(structure, chain, seqid);
    if (!potential_site.has_value()) {
        throw std::runtime_error("Site could not be found");
    }
    Sails::Glycosites glycosites = {potential_site.value()};
    return run_em_cycle(glycosites, structure, grid, resolution, cycles, resource_dir, false, verbose);
}



//SNFG FUNCTIONS

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


std::map<std::string, std::string> find_wurcs(gemmi::Structure& structure, std::string& chain, int seqid, std::string& resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    Sails::Topology topology = {&structure, residue_database};

    std::optional<Sails::Glycosite> potential_glycosite = Sails::find_site(structure, chain, seqid);
    if (!potential_glycosite.has_value()) throw std::runtime_error("Could not find specified site");
    Sails::Glycosite glycosite = potential_glycosite.value();
    Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
    std::string generated_wurcs =  Sails::WURCS::generate_wurcs(&glycan, residue_database);
    std::string key = Sails::Utils::format_residue_from_site(glycosite, &structure);
    std::map<std::string, std::string> wurcs_map = {{key, generated_wurcs}};

    return wurcs_map;
}


std::map<std::string, std::string> find_all_wurcs(gemmi::Structure& structure, std::string& resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    Sails::Topology topology = {&structure, residue_database};

    std::map<std::string, std::string> wurcs_map;

    Sails::Glycosites n_glycosites = Sails::find_n_glycosylation_sites(structure);
    for (auto& site: n_glycosites) {
        Sails::Glycan glycan = topology.find_glycan_topology(site);
        if (glycan.empty()) continue;
        std::string generated_wurcs =  Sails::WURCS::generate_wurcs(&glycan, residue_database);
        std::string key = Sails::Utils::format_residue_from_site(site, &structure);
        wurcs_map[key] = generated_wurcs;
    }

    Sails::Glycosites c_glycosites = Sails::find_c_glycosylation_sites(structure);
    for (auto& site: c_glycosites) {
        Sails::Glycan glycan = topology.find_glycan_topology(site);
        if (glycan.empty()) continue;
        std::string generated_wurcs =  Sails::WURCS::generate_wurcs(&glycan, residue_database);
        std::string key = Sails::Utils::format_residue_from_site(site, &structure);
        wurcs_map[key] = generated_wurcs;
    }

    //
    // if (!potential_glycosite.has_value()) throw std::runtime_error("Could not find specified site");
    // Sails::Glycosite glycosite = potential_glycosite.value();
    // Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
    // std::string generated_wurcs =  Sails::WURCS::generate_wurcs(&glycan, residue_database);

    return wurcs_map;
}

gemmi::Structure model_wurcs(gemmi::Structure& structure, std::string& wurcs, std::string& chain, int seqid, std::string& resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    std::optional<Sails::Glycosite> potential_glycosite = Sails::find_site(structure, chain, seqid);
    if (!potential_glycosite.has_value()) throw std::runtime_error("Could not find specified site");
    Sails::Glycosite glycosite = potential_glycosite.value();

    Sails::Model model = {&structure, linkage_database, residue_database};
    Sails::PseudoGlycan pseudo_glycan  = Sails::WURCS::generate_pseudo_glycan(wurcs, &structure, glycosite, linkage_database, residue_database);
    model.create_pseudo_glycan(pseudo_glycan);

    return structure;
}


gemmi::Structure morph(gemmi::Structure& structure, std::string& wurcs, std::string& chain, int seqid, std::string& resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    std::optional<Sails::Glycosite> potential_glycosite = Sails::find_site(structure, chain, seqid);
    if (!potential_glycosite.has_value()) throw std::runtime_error("Could not find specified site");
    Sails::Glycosite glycosite = potential_glycosite.value();

    Sails::Topology topology = {&structure, residue_database};

    Sails::Glycan glycan = topology.find_glycan_topology(glycosite);

    Sails::PseudoGlycan pseudo_glycan  = Sails::WURCS::generate_pseudo_glycan(wurcs, &structure, glycosite, linkage_database, residue_database);
    Sails::Morph morpher = {&structure};
    morpher.transform(glycan, pseudo_glycan);

    return structure;
}


Sails::Output validate(gemmi::Structure& structure, Sails::MTZ &sails_mtz, bool remove, float threshold, std::string& resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    gemmi::Mtz mtz = form_gemmi_mtz(sails_mtz);
    check_spacegroup(&mtz, &structure); // check to ensure the MTZ has a spacegroup

    auto density = Sails::XtalDensity(mtz);
    density.load_map_coefficients();

    std::map<Sails::Glycosite, double> rsccs = Sails::Score::calculate_rsccs(&density, &structure, residue_database);

    std::vector<Sails::Glycosite> to_remove = {};
    std::vector<Sails::TelemetryFormat> log = {};

    for (auto& [site, rscc]: rsccs) {
        std::string residue_key = Sails::Utils::format_residue_from_site(site, &structure);
        log.emplace_back(residue_key, rscc);
        if (rscc > threshold) {
                continue;
        }
        to_remove.emplace_back(site);
    }

    if (remove) {
        Sails::Topology topology = {&structure, residue_database};

        std::set<Sails::Glycosite> removal_set = {to_remove.begin(), to_remove.end()};

        for (auto &site: to_remove) {
            auto glycan = topology.find_glycan_topology(site);
            std::vector<Sails::Sugar*> downstream_sugars = glycan.get_downstream_sugars(site);
            for (auto& downstream_sugar: downstream_sugars) {
                if (std::find(removal_set.begin(), removal_set.end(), downstream_sugar->site) != removal_set.end()) continue;
                downstream_sugar->site.atom_idx = 0; // remove atom site from site to allow sorting
                removal_set.insert(downstream_sugar->site);
            }
        }

        std::vector<Sails::Glycosite> removal_list = {removal_set.begin(), removal_set.end()};

        std::sort(removal_list.begin(), removal_list.end(), [](const Sails::Glycosite& a, const Sails::Glycosite& b) {
            return !(a < b);
        });

        for (auto &site: removal_list) {
            const auto residue_ptr = &structure.models[site.model_idx].chains[site.chain_idx].residues;
            residue_ptr->erase(residue_ptr->begin() + site.residue_idx);
        }
    }


    std::string log_string = Sails::Telemetry::format_log(log, false, "").value();
    return {
        structure,
        log_string
    };
}

Sails::Output validate(gemmi::Structure& structure, gemmi::Grid<>& grid, float resolution, bool remove, float threshold, bool use_q, std::string& resource_dir) {
    std::string data_file = resource_dir + "/data.json";
    Sails::JSONLoader loader = {data_file};
    Sails::ResidueDatabase residue_database = loader.load_residue_database();
    Sails::LinkageDatabase linkage_database = loader.load_linkage_database();

    auto density = Sails::EMDensity(grid, resolution);


    std::map<Sails::Glycosite, double> rsccs = Sails::Score::calculate_rsccs(&density, &structure, residue_database);
    std::map<Sails::Glycosite, double> qscores = Sails::Score::calculate_qscores(&density, &structure, residue_database);
    std::map<Sails::Glycosite, double> scores = use_q ? qscores : rsccs;

    // equation from https://doi.org/10.1107/S2059798325005923
    double q_score_threshold = -0.0016*pow(resolution,2) + 0.0434*pow(resolution,2)-0.3956*resolution + 1.3366;

    double applied_threshold = use_q ? q_score_threshold : threshold ;

    if (remove) {
        std::cout << "Enforcing score limit of " << applied_threshold << std::endl;
    }
    std::vector<Sails::Glycosite> to_remove = {};
    std::vector<Sails::TelemetryFormat> log = {};

    for (auto& [site, score]: scores) {
        std::string residue_key = Sails::Utils::format_residue_from_site(site, &structure);
        log.emplace_back(residue_key, rsccs.at(site), qscores.at(site));
        if (score > applied_threshold) {
                continue;
        }
        to_remove.emplace_back(site);
    }

    if (remove) {
        Sails::Topology topology = {&structure, residue_database};

        std::set<Sails::Glycosite> removal_set = {to_remove.begin(), to_remove.end()};

        for (auto &site: to_remove) {
            auto glycan = topology.find_glycan_topology(site);
            std::vector<Sails::Sugar*> downstream_sugars = glycan.get_downstream_sugars(site);
            for (auto& downstream_sugar: downstream_sugars) {
                if (removal_set.count(downstream_sugar->site) > 0) continue;
                downstream_sugar->site.atom_idx = 0; // remove atom site from site to allow sorting
                removal_set.insert(downstream_sugar->site);
            }
        }

        std::vector<Sails::Glycosite> removal_list = {removal_set.begin(), removal_set.end()};

        std::sort(removal_list.begin(), removal_list.end(), [](const Sails::Glycosite& a, const Sails::Glycosite& b) {
            return !(a < b);
        });

        for (auto &site: removal_list) {
            const auto residue_ptr = &structure.models[site.model_idx].chains[site.chain_idx].residues;
            residue_ptr->erase(residue_ptr->begin() + site.residue_idx);
        }
    }


    std::string log_string = Sails::Telemetry::format_log(log, false, "").value();
    return {
        structure,
        log_string
    };
}


// gemmi::Structure wurcs(gemmi::Structure& structure, std::string chain, int seqid, std::string& resource_dir) {
//     std::string data_file = resource_dir + "/data.json";
//     Sails::JSONLoader loader = {data_file};
//     Sails::ResidueDatabase residue_database = loader.load_residue_database();
//     Sails::LinkageDatabase linkage_database = loader.load_linkage_database();
//
//     // Sails::Topology topology = {&structure, residue_database};
//     //
//     std::optional<Sails::Glycosite> potential_glycosite = Sails::find_site(structure, chain, seqid);
//     if (!potential_glycosite.has_value()) throw std::runtime_error("Could not find specified site");
//     Sails::Glycosite glycosite = potential_glycosite.value();
//     //
//     // Sails::Glycan glycan = topology.find_glycan_topology(glycosite);
//     //
//     // std::string generated_wurcs =  Sails::WURCS::generate_wurcs(&glycan, residue_database);
//
//     // std::cout << "GENERATED WURCS: " << generated_wurcs << std::endl;
//     // std::string wurcs = "WURCS=2.0/3,6,5/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2/a4-b1_b4-c1";
//     std::string wurcs = "WURCS=2.0/3,6,5/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3-3/a4-b1_b4-c1_c3-d1_d2-e1_e2-f1";
//     Sails::WURCS::generate_pseudo_glycan(wurcs, &structure, glycosite, linkage_database, residue_database);
//
//     return structure;
// }


void test() {
//    const std::string path = "testing/test_data/4ax7/4AX7_deglycosylated.cif";
//    const std::string mtz_path = "testing/test_data/4ax7/4AX7.mtz";
//    gemmi::Mtz mtz = gemmi::read_mtz_file(mtz_path);
//    auto smtz = Sails::form_sails_mtz(mtz, "FP", "SIGFP");
//    gemmi::Structure structure = gemmi::read_structure_file(path);
//
//    std::string data_file = "package/src/sails/data/data.json";
//    Sails::JSONLoader loader = {data_file};
//    Sails::ResidueDatabase residue_database = loader.load_residue_database();
//
//    Sails::Density density = Sails::Density(mtz);
//    density.load_hkl("FP", "SIGFP");
//    density.recalculate_map(structure);


//    auto o = find_o_mannosylation_sites(structure, sa_map);
//    std::string a = "package/src/sails/data";
//    auto output = run_cycle(o, structure, smtz, 1, a, true, true);
//    Sails::Utils::save_structure_to_file(output.structure, "o-mannose-strict.cif");
//    std::cout << output.log << std::endl;
}

// testbed
int main() {
    const std::string path = "testing/em/7lze/7lze.cif";
    const std::string mtz_path = "testing/em/7lze/emd_23605.map";

    gemmi::Structure structure = gemmi::read_structure_file(path);
    gemmi::Ccp4<float> map;
    map.read_ccp4_file(mtz_path);

    std::string data_file = "package/src/sails/data/";
    auto glycosites = Sails::find_n_glycosylation_sites(structure);

    // run_em_cycle(glycosites, structure, map.grid, 1, data_file,  false, true);
}
