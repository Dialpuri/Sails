//
// Created by jordan on 05/11/24.
//

#include "../include/sails-wurcs.h"

std::string Sails::WURCS::generate_wurcs(Sails::Glycan *glycan, Sails::ResidueDatabase& residue_database) {

    std::stringstream wurcs;

    wurcs << get_wurcs_header() << get_section_delimiter();
    wurcs << get_unit_count(glycan) << get_section_delimiter();
    wurcs << get_unique_residue_list(glycan, residue_database) << get_section_delimiter();

    return wurcs.str();
}

Sails::Glycan Sails::WURCS::generate_psuedo_glycan(const std::string &wurcs, gemmi::Structure *structure) {
    return Sails::Glycan();
}

std::string Sails::WURCS::get_unit_count(Sails::Glycan *glycan) {
    std::stringstream unit_count;

    unit_count << glycan->unique_sugar_count() << ",";
    unit_count << glycan->sugar_count() << ",";
    unit_count << glycan->linkage_count();
    return unit_count.str();
}

std::string Sails::WURCS::get_unique_residue_list(Sails::Glycan *glycan, Sails::ResidueDatabase &residue_database) {

    std::vector<std::string> unique_sugars = glycan->get_unique_sugar_names();

    std::stringstream unique_residues;
    for (auto& sugar: unique_sugars) {
        std::cout << sugar << std::endl;
        if (residue_database.find(sugar) == residue_database.end()) {
            unique_residues << "[" << sugar << "]";
            continue;
        }
        std::optional<std::string> wurcs = residue_database[sugar].wurcs_code;
        if (!wurcs.has_value()) {continue;} // an ASN, TRP, etc. would not have a WURCS code

        unique_residues << "[" << wurcs.value() << "]";
    }

    return unique_residues.str();
}
