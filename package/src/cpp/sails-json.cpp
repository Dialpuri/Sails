//
// Created by jordan on 05/07/24.
//

#include "../include/sails-json.h"


Sails::JSONLoader::JSONLoader(const std::string& path, Sails::ResidueDatabase& database) {
    this->init(path, database);
}

void Sails::JSONLoader::init(const std::string &path, Sails::ResidueDatabase& database) {
    if (!Sails::Utils::check_file_exists(path)) {throw std::runtime_error("Data file could not be found");}

    simdjson::ondemand::parser parser;
    auto json = simdjson::padded_string::load(path);
    simdjson::ondemand::document doc = parser.iterate(json);

    // load residues
    auto residues = doc["residues"];
    for (auto value: residues) {
        std::string name = std::string(value["name"].get_string().value());

        auto donor_sets = extract_atom_set(value, "donor_sets");
        auto acceptors_sets = extract_atom_set(value, "acceptor_sets");

        ResidueData data = {acceptors_sets, donor_sets};
        database.insert({name,  data});
    }

}


std::vector<Sails::AtomSet> Sails::JSONLoader::extract_atom_set(simdjson::simdjson_result<simdjson::ondemand::value> &value, const char* key) const {
    simdjson::simdjson_result<simdjson::ondemand::array> sets = value[key].get_array();

    std::vector<Sails::AtomSet> atom_sets;
    for (auto set: sets) {
        std::string atom_1 = std::string(set["atom_1"].get_string().value());
        std::string atom_2 = std::string(set["atom_2"].get_string().value());
        std::string atom_3 = std::string(set["atom_3"].get_string().value());
        int identifier = set["identifier"].get_int64();

        Sails::AtomSet atom_set = {atom_1, atom_2, atom_3, identifier};
        atom_sets.emplace_back(atom_set);
    }

    return atom_sets;
}
