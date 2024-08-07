//
// Created by Jordan Dialpuri on 31/07/2024.
//

#include "../include/sails-telemetry.h"

#include <src/include/sails-json.h>


void Sails::Telemetry::save_state(int cycle) {
    if (cycle == 1) {
        states[cycle] = sites;
        return;
    }
    std::set<Glycosite> difference;
    std::set_difference(sites.begin(), sites.end(), states[cycle - 1].begin(), states[cycle - 1].end(),
                        std::inserter(difference, difference.begin()));
    states[cycle] = difference;
}

void Sails::Telemetry::save_snfg(int cycle, std::string& key, std::string &snfg) {
    snfgs[cycle][key] = snfg;
}

void Sails::Telemetry::format_log(gemmi::Structure *structure) {
    for (const auto &[cycle, sites]: states) {
        std::cout << "Cycle " << cycle << "\n";
        for (const auto &site: sites) {
            std::cout << "Added " << Utils::format_residue_from_site(site, structure) << "\n";
        }
        std::cout << std::endl;
    }
}

Sails::TelemetryLog Sails::Telemetry::calculate_log(gemmi::Structure *structure, Density *density) {
    TelemetryLog log;
    for (const auto &[cycle, sites]: states) {
        for (const auto &site: sites) {
            gemmi::Residue residue = Utils::get_residue_from_glycosite(site, structure);
            const double rscc_score = density->score_residue(residue, rscc);
            const double rsr_score = density->score_residue(residue, rsr);
            const double dds_score = density->score_residue(residue, dds);
            log[cycle].emplace_back(
                Utils::format_residue_from_site(site, structure),
                rscc_score,
                rsr_score,
                dds_score);
        }
    }
    return log;
}


std::optional<std::string> Sails::Telemetry::format_log(gemmi::Structure *structure, Density *density, bool write) {
    TelemetryLog log = calculate_log(structure, density);

    JSONWriter writer;
    if (write) {
        std::ofstream stream(m_filepath);
        writer.write_json_file(log, stream);
        stream.close();
    } else {
        std::stringstream stream;
        writer.write_json_file(log, stream);
        return stream.str();
    }
    return std::nullopt;

    // for (const auto &[cycle, entries]: log) {
    //     std::cout << "Cycle " << cycle << "\n";
    //     for (const auto &entry: entries) {
    //         std::cout << "Added " << entry.residue_id << "\t" << entry.rscc_score << "\t" << entry.rsr_score << "\t"
    //         << entry.dds_score << "\n";
    //     }
    //     std::cout << std::endl;
    // }
}
