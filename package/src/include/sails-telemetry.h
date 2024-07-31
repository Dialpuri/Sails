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
    /**
     * @class Telemetry
     *
     * @brief Handles telemetry for Sails addition and deletions.
     *
     * This class provides functionality to store and manipulate telemetry data
     * related to Glycosite objects. It supports operations such as adding and
     * removing Glycosite objects, retrieving the size of the telemetry data,
     * saving the state of the telemetry, and formatting the log with a given
     * gemmi::Structure object.
     */
    struct Telemetry {
        explicit Telemetry(const std::string& filepath): filepath_(filepath) {}

        /**
         * @brief Overloaded operator<< for inserting Glycosite objects into a set.
         *
         * This function takes a Glycosite object and inserts it into a set of Glycosite objects.
         * The set is used to store telemetry data related to Glycosite objects.
         * This operator is typically used in conjunction with std::cout to output Glycosite objects.
         *
         * @param site The Glycosite object to be inserted.
         */
        void operator<<(const Glycosite& site) {
            sites.insert(site);
        };

        /**
         * @brief Overloaded operator>> to remove a Glycosite from the sites container.
         *
         * This operator removes the provided Glycosite object from the sites container.
         *
         * @param site The Glycosite object to be removed.
         */
        void operator>>(const Glycosite& site) {
            sites.erase(site);
        }

        /**
         * @brief Returns the size of the telemetry data.
         *
         * This method returns the number of elements in the "sites" container
         * which represents the telemetry data. It provides a way to determine
         * the current size of the telemetry data stored in the container.
         *
         * @return The number of elements in the telemetry data.
         */
        [[nodiscard]] size_t size() const {return sites.size();}

        /**
         * @brief Saves the state of the telemetry at the given cycle.
         *
         * This method saves the state of the telemetry at the specified cycle
         * by calculating the difference between the current set of Glycosite objects
         * and the set of Glycosite objects at the previous cycle.
         * If the given cycle is greater than 1, it calculates the difference using
         * std::set_difference() and stores the result in the states map at the given cycle.
         * If the given cycle is 1 or less, it simply stores the current set of Glycosite objects
         * in the states map at the given cycle.
         *
         * @param cycle The cycle at which to save the state of the telemetry.
         */
        void save_state(int cycle);

        /**
         * @brief Formats the log with a given gemmi::Structure object.
         *
         * This method prints the formatted log of the telemetry data,
         * including the cycles and added Glycosite objects in each cycle.
         * The log is printed to the standard output.
         *
         * @param structure A pointer to the gemmi::Structure object used for formatting the log.
         *                  This object is used to retrieve information about Glycosite objects.
         */
        void format_log(gemmi::Structure* structure);

    private:
        std::set<Glycosite> sites{};
        std::map<int, std::set<Glycosite>> states{};
        const std::string &filepath_;
    };

}
#endif //SAILS_TELEMETRY_H
