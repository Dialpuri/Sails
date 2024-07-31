//
// Created by Jordan Dialpuri on 31/07/2024.
//

#ifndef SAILS_TELEMETRY_H
#define SAILS_TELEMETRY_H

#include <iostream>
#include <set>
#include "sails-model.h"
#include "sails-utils.h"

namespace Sails {

    struct Telemetry {
        explicit Telemetry(const std::string& filepath): filepath_(filepath) {}

        void operator<<(const Glycosite& site) {
            sites.insert(site);
        };

        void operator>>(const Glycosite& site) {
            sites.erase(site);
        }

        [[nodiscard]] size_t size() const {return sites.size();}

        void save_state(int cycle);

        void format_log(gemmi::Structure* structure);

    private:
        std::set<Glycosite> sites{};
        std::map<int, std::set<Glycosite>> states{};
        const std::string &filepath_;
    };

}
#endif //SAILS_TELEMETRY_H
