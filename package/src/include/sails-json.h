//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_JSON_H
#define SAILS_SAILS_JSON_H

#include "../third-party/simdjson.h"
#include "sails-model.h"

#include <filesystem>

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
        JSONLoader(const std::string& path, Sails::ResidueDatabase& database);
        void init(const std::string& path, Sails::ResidueDatabase& database);

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
        std::vector<Sails::AtomSet> extract_atom_set(simdjson::simdjson_result<simdjson::ondemand::value> &value, const char *key) const;


    }; // class JSONLoader

} // namespace Sails

#endif //SAILS_SAILS_JSON_H
