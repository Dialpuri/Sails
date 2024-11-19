//
// Created by jordan on 05/11/24.
//

#include "../include/sails-wurcs.h"

std::string Sails::WURCS::generate_wurcs(Sails::Glycan *glycan, Sails::ResidueDatabase& residue_database) {

    std::stringstream wurcs;

    wurcs << get_wurcs_header() << get_section_delimiter();
    wurcs << get_unit_count(glycan) << get_section_delimiter();
    wurcs << get_unique_residue_list(glycan, residue_database) << get_section_delimiter();
    wurcs << get_residue_order(glycan, residue_database) << get_section_delimiter();
    wurcs << get_link_list(glycan, residue_database);
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

std::string Sails::WURCS::get_residue_order(Sails::Glycan *glycan, Sails::ResidueDatabase &residue_database) {

    std::vector<int> order = calculate_residue_order(glycan, residue_database);

    std::stringstream residue_order;
    for (int i = 0; i < order.size(); i++) {
        residue_order << order[i];

        if (i != order.size() - 1) {
            residue_order << "-";
        }
    }

    return residue_order.str();
}

std::string Sails::WURCS::get_link_list(Sails::Glycan *glycan, Sails::ResidueDatabase &residue_database) {

    std::vector<int> order = calculate_residue_order(glycan, residue_database);
    auto sugar_sites = glycan->get_sugar_site_dfs_order();
    auto sugar_names = glycan->get_sugar_name_order();

    auto sugars = glycan->get_sugars();
    auto adjacency_list = glycan->get_adjacency_list();

    char alphabet[] = "abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::map<Glycosite, char> site_map;
    int alphabet_count = 0;

    for (int i = 0; i < sugar_sites.size(); i++) {
        if (residue_database.find(sugar_names[i]) == residue_database.end()) {
            continue;
        }
        std::optional<std::string> wurcs = residue_database[sugar_names[i]].wurcs_code;
        if (!wurcs.has_value()) {continue;}

        site_map[sugar_sites[i]] = alphabet[alphabet_count];
        alphabet_count++;
    }

    std::vector<std::string> linkage_order;

    for (auto& site: sugar_sites) {
        Sugar* sugar_ptr = sugars->at(site).get();
        if (adjacency_list->find(sugar_ptr) == adjacency_list->end()) {
            continue;
        }

        std::set<Sugar*> linked_sugars = adjacency_list->at(sugar_ptr);

        std::map<int, Linkage*> linkages;
        for (auto& linked_sugar: linked_sugars) {
            Linkage* linkage = sugar_ptr->find_linkage(sugar_ptr, linked_sugar);

            if (site_map.find(linkage->donor_sugar->site) == site_map.end()) {
                continue;
            }

            if (site_map.find(linkage->acceptor_sugar->site) == site_map.end()) {
                continue;
            }

            int linkage_donor_number = linkage->donor_atom.back() - '0';
            linkages[linkage_donor_number] = linkage;
        }

        for (auto& [number, linkage]: linkages) {
            std::stringstream linkage_string;
            linkage_string << site_map.at(linkage->donor_sugar->site) << number;
            linkage_string << "-";
            linkage_string << site_map.at(linkage->acceptor_sugar->site) << linkage->acceptor_atom.back();
            linkage_order.emplace_back(linkage_string.str());
        }
    }

    std::stringstream linkage_stream;
    for (int i = 0; i < linkage_order.size(); i++) {
        linkage_stream << linkage_order[i];
        if (i < linkage_order.size()-1) {
            linkage_stream << "_";
        }
    }

    return linkage_stream.str();
}

std::vector<int>
Sails::WURCS::calculate_residue_order(Sails::Glycan *glycan, Sails::ResidueDatabase &residue_database) {
    std::vector<std::string> sugar_order = glycan->get_sugar_name_dfs_order();
    std::vector<std::string> unique_sugars = glycan->get_unique_sugar_names();
    std::map<std::string, int> sugar_indices;

    for (int i = 0; i < unique_sugars.size(); i++) {
        std::string unique_sugar = unique_sugars[i];
        std::optional<std::string> wurcs = residue_database[unique_sugar].wurcs_code;
        if (!wurcs.has_value()) {continue;}
        sugar_indices[unique_sugars[i]] = i;
    }

    std::vector<int> order ;
    for (const auto& name: sugar_order) {
        if (sugar_indices.find(name) == sugar_indices.end()) {
            continue;
        }
        int index = sugar_indices[name];
        order.emplace_back(index);
    }
    return order;
}
