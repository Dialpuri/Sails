//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_JSON_H
#define SAILS_SAILS_JSON_H

#include "../third-party/simdjson.h"
#include "sails-model.h"

#include <filesystem>

#include "sails-telemetry.h"

namespace Sails {

    /**
     * @class JSONLoader
     * @brief The JSONLoader class is responsible for loading and parsing JSON files.
     *
     * This class provides methods for loading a JSON data file and extracting data from it.
     * It requires a path to the file and ResidueDatabase object to fill.
     */
    class JSONLoader {

    public:
        JSONLoader() = default;

        JSONLoader(const std::string &path);

        void init(const std::string &path);

        /**
         * @brief Loads the residue database from a JSON file.
         *
         * This method loads a JSON file containing residue data and parses it to create a ResidueDatabase object.
         * It uses the simdjson library for efficient parsing. The path to the JSON file must be specified in the 'm_path' member variable.
         *
         * @return The ResidueDatabase object containing the loaded residue data.
         */
        ResidueDatabase load_residue_database();

        /**
        * @brief Loads the linkage database from a JSON file.
        *
        * This method loads a JSON file containing linkage data and parses it to create a LinkageDatabase object.
        * It uses the simdjson library for efficient parsing. The path to the JSON file must be specified in the 'm_path' member variable.
        *
        * @return The LinkageDatabase object containing the loaded residue data.
        */
        LinkageDatabase load_linkage_database();

    private:
        /**
         * @brief Extracts a vector of AtomSet objects from a JSON value using the specified key.
         *
         * This function takes a simdjson value and extracts an array of objects from it using the specified key.
         * Each object in the array represents an AtomSet, which contains information about three atoms and an identifier.
         * The extracted AtomSets are returned as a vector.
         *
         * @param value The simjson value to extract data from.
         * @param key The key to extract the array of AtomSet objects.
         * @return std::vector<Sails::AtomSet> The vector of extracted AtomSet objects.
         */
        static std::vector<AtomSet>
        extract_atom_set(simdjson::simdjson_result<simdjson::ondemand::value> &value, const char *key);

        /**
         * @brief Extracts angles from a JSON object.
         *
         * This method extracts angles from a JSON object using keys "alpha", "beta", and "gamma".
         *
         * @param value A simdjson_result<simdjson::ondemand::value> reference to extract angles from.
         *
         * @return An AngleSet object containing the extracted angles.
         */
        static AngleSet extract_angles(simdjson::simdjson_result<simdjson::ondemand::value> &value);


        /**
         * @brief Extracts the torsion angles from a JSON value.
         *
         * This method extracts angles from a JSON object using the keys "phiMean", "phiStdDev", "psiMean",
         * "psiStdDev", "omegaMean", and "omegaStdDev", each associated with a numerical value.
         *
         * @param value A simdjson_result<simdjson::ondemand::value> reference to extract torsions from..
         * @return A TorsionSet object containing the extracted torsion angles.
         */
        static TorsionSet extract_torsions(simdjson::simdjson_result<simdjson::ondemand::value> &value);


    private:
        std::string m_path;
        simdjson::ondemand::document m_doc;

    }; // class JSONLoader


    /**
     * @class JSONWriter
     * @brief The JSONWriter class is responsible for writing telemetry log data to a JSON file.
     *
     * This class provides a method for writing telemetry log data to a JSON file. It requires a
     * filename to write the data to. The class uses the TelemetryLog object to access the log data
     * and construct the JSON structure.
     */
    class JSONWriter {
    public:
        JSONWriter(const std::string& filename): m_filename(filename) {}

        /**
         * @brief Writes the telemetry log data to a JSON file.
         *
         * This method takes a TelemetryLog object and writes its data to a JSON file.
         * The JSON file includes information such as the date, cycles, and entries.
         *
         * @param log The TelemetryLog object containing the telemetry data to write.
         *
         * @note The log object should be populated with data prior to calling this method.
         * @note The filename of the JSON file to write is specified in the constructor of the JSONWriter class.
         *
         * @see TelemetryLog
         * @see JSONWriter
         */
        void write_json_file(TelemetryLog& log);

    private:
        std::string m_filename;
    };

} // namespace Sails

#endif //SAILS_SAILS_JSON_H
