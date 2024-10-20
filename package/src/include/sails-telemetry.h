//
// Created by Jordan Dialpuri on 31/07/2024.
//

#ifndef SAILS_TELEMETRY_H
#define SAILS_TELEMETRY_H

#include <iostream>
#include <set>
#include <optional>

#include "sails-model.h"
#include "sails-utils.h"
#include "density/sails-density.h"

namespace Sails {

    struct TelemetryFormat {
        TelemetryFormat() = default;

        TelemetryFormat(const std::string &residue_id, double rscc_score, double rsr_score, double dds_score)
            : residue_id(residue_id),
              rscc_score(rscc_score),
              rsr_score(rsr_score),
              dds_score(dds_score) {
        }

        std::string residue_id;
        double rscc_score;
        double rsr_score;
        double dds_score;
    };
    typedef std::map<int, std::vector<TelemetryFormat>> TelemetryLog;


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
        explicit Telemetry(const std::string& filepath): m_filepath(filepath) {}

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
         *
         * @brief Overloaded operator<< for appending multiple Glycosite objects to the stream.
         *
         * This function appends multiple Glycosite objects from the given vector to the stream.
         * It iterates through the vector and calls the operator<<(const Glycosite& site) for each Glycosite object.
         *
         * @param sites A vector of Glycosite objects to be appended to the stream.
         */
        void operator<<(const std::set<Glycosite>& sites) {
            for (const auto& site: sites) {
                this->operator<<(site);
            }
        }


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
         *
         * @brief Overloaded operator>> for appending multiple Glycosite objects to the stream.
         *
         * This function appends multiple Glycosite objects from the given vector to the stream.
         * It iterates through the vector and calls the operator>>(const Glycosite& site) for each Glycosite object.
         *
         * @param sites A vector of Glycosite objects to be appended to the stream.
         */
        void operator>>(const std::set<Glycosite>& sites) {
            for (const auto& site: sites) {
                this->operator>>(site);
            }
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
         * If the given cycle is not 1 it calculates the difference using
         * std::set_difference() and stores the result in the states map at the given cycle.
         * If the given cycle is 1, it simply stores the current set of Glycosite objects
         * in the states map at the given cycle.
         *
         * @param cycle The cycle at which to save the state of the telemetry.
         */
        void save_state(int cycle);


        /**
         * @brief Saves the SNFG (Symbol Nomenclature for Glycans) representation of a glycan structure.
         *
         * This method saves the SNFG representation of a glycan structure for a specific cycle and glycosite.
         * The SNFG representation is stored as a string in the provided `snfg` parameter.
         *
         * @param cycle The cycle number of the glycan structure.
         * @param key The key corresponding to the glycosite.
         * @param snfg The string to store the SNFG representation of the glycan structure.
         */
        void save_snfg(int cycle, std::string& key, std::string &snfg);

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


        /**
         * @brief Calculates the telemetry log for Sails.
         *
         * This method takes a gemmi::Structure pointer and a Density pointer and calculates
         * the telemetry log for Sails. It iterates through the states and sites and calculates
         * the rscc, rsr, and dds scores for each site by calling the density->score_residue method.
         * It then formats the residue and adds it along with the scores to the telemetry log.
         * Finally, it returns the telemetry log.
         *
         * @param structure A pointer to the gemmi::Structure object.
         * @param density A pointer to the Density object.
         * @return The calculated telemetry log for Sails.
         */
        TelemetryLog calculate_log(gemmi::Structure *structure, Density *density);

        /**
         * @brief Formats the log with the given gemmi::Structure object.
         *
         * This method formats the log by printing information about the cycles and
         * added Glycosite objects stored in the telemetry. It uses the provided
         * gemmi::Structure object to retrieve residue and score information for each
         * Glycosite.
         *
         * @param structure A pointer to the gemmi::Structure object.
         * @param density A pointer to the Density object.
         * @param write
         *
         * @note The method assumes that the telemetry data has already been stored
         * using the add_glycosite() method.
         *
         * @see Sails::Telemetry::add_glycosite()
         * @see Utils::format_residue_from_site()
         * @see Utils::get_residue_from_glycosite()
         * @see Density::score_residue()
         */
        std::optional<std::string> format_log(gemmi::Structure *structure, Density* density, bool write);

        typedef std::map<int, std::map<std::string, std::string>> SNFGCycleData;
        /**
         * @brief Get the SNFGs.
         *
         * This method returns the SNFGs (Symbol Nomenclature for Glycans) as a map of maps.
         * The outer map uses the cycle number as keys.
         * The inner map uses the base residue as keys and the SNFGs are the values
         *
         * @return A map of maps representing the SNFGs. The outer map uses integers as keys
         *         and the inner map uses strings as keys.
         */
        [[nodiscard]] SNFGCycleData get_snfgs() const {return snfgs;};


    private:
        std::set<Glycosite> sites{};
        std::map<int, std::set<Glycosite>> states{};
        SNFGCycleData snfgs;
        const std::string &m_filepath;
    };

}
#endif //SAILS_TELEMETRY_H
