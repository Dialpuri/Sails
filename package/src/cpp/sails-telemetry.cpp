//
// Created by Jordan Dialpuri on 31/07/2024.
//

#include "../include/sails-telemetry.h"

void Sails::Telemetry::save_state(int cycle) {
    if (cycle > 1) {
        std::set<Glycosite> difference;
        std::set_difference(sites.begin(), sites.end(), states[cycle-1].begin(), states[cycle-1].end(),
                            std::inserter(difference, difference.begin()));
        states[cycle] = difference;
    } else {
        states[cycle] = sites;
    }
}

void Sails::Telemetry::format_log(gemmi::Structure *structure) {
    for (const auto& [cycle, sites]: states) {
        std::cout << "Cycle " << cycle << "\n";
        for (const auto& site: sites) {
            std::cout << "Added " << Utils::format_residue_from_site(site, structure) << "\n";
        }
        std::cout << std::endl;
    }
}
