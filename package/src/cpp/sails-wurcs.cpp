//
// Created by jordan on 05/11/24.
//

#include "../include/sails-wurcs.h"
#include "../include/sails-linkage.h"

// WURCS WRITING

std::string Sails::WURCS::generate_wurcs(Sails::Glycan *glycan, Sails::ResidueDatabase &residue_database) {
    std::stringstream wurcs;

    wurcs << get_wurcs_header() << get_section_delimiter();
    wurcs << get_unit_count(glycan) << get_section_delimiter();
    wurcs << get_unique_residue_list(glycan, residue_database) << get_section_delimiter();
    wurcs << get_residue_order(glycan, residue_database) << get_section_delimiter();
    wurcs << get_link_list(glycan, residue_database);
    return wurcs.str();
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
    for (auto &sugar: unique_sugars) {
        if (residue_database.find(sugar) == residue_database.end()) {
            unique_residues << "[" << sugar << "]";
            continue;
        }
        std::optional<std::string> wurcs = residue_database[sugar].wurcs_code;
        if (!wurcs.has_value()) { continue; } // an ASN, TRP, etc. would not have a WURCS code

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
        if (!wurcs.has_value()) { continue; }

        site_map[sugar_sites[i]] = alphabet[alphabet_count];
        alphabet_count++;
    }

    std::vector<std::string> linkage_order;

    for (auto &site: sugar_sites) {
        Sugar *sugar_ptr = sugars->at(site).get();
        if (adjacency_list->find(sugar_ptr) == adjacency_list->end()) {
            continue;
        }

        std::set<Sugar *> linked_sugars = adjacency_list->at(sugar_ptr);

        std::map<int, Linkage *> linkages;
        for (auto &linked_sugar: linked_sugars) {
            Linkage *linkage = sugar_ptr->find_linkage(sugar_ptr, linked_sugar);

            if (site_map.find(linkage->donor_sugar->site) == site_map.end()) {
                continue;
            }

            if (site_map.find(linkage->acceptor_sugar->site) == site_map.end()) {
                continue;
            }

            int linkage_donor_number = linkage->donor_atom.back() - '0';
            linkages[linkage_donor_number] = linkage;
        }

        for (auto &[number, linkage]: linkages) {
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
        if (i < linkage_order.size() - 1) {
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
        if (!wurcs.has_value()) { continue; }
        sugar_indices[unique_sugars[i]] = i;
    }

    std::vector<int> order;
    for (const auto &name: sugar_order) {
        if (sugar_indices.find(name) == sugar_indices.end()) {
            continue;
        }
        int index = sugar_indices[name];
        order.emplace_back(index);
    }
    return order;
}

// WURCS PARSING


std::vector<std::string> Sails::WURCS::extract_wurcs_unique_residues(const std::string &wurcs) {
    std::regex regex(R"(\[\S*\])");
    std::vector<std::string> results;
    auto matches_begin = std::sregex_iterator(wurcs.begin(), wurcs.end(), regex);
    auto matches_end = std::sregex_iterator();

    for (auto it = matches_begin; it != matches_end; ++it) {
        results.push_back(it->str());
    }

    if (results.empty() || results.size() > 1) { throw std::runtime_error("Error in WURCS Unique Residue List"); }

    std::vector<std::string> order = Utils::split(results[0], ']');
    for (auto &i: order) {
        i.erase(0, 1);
    }

    return order;
}

std::vector<int> Sails::WURCS::extract_wurcs_residue_order(const std::string &wurcs) {
    std::regex regex("/[0-9-]*/");
    std::vector<std::string> results;
    auto matches_begin = std::sregex_iterator(wurcs.begin(), wurcs.end(), regex);
    auto matches_end = std::sregex_iterator();

    for (auto it = matches_begin; it != matches_end; ++it) {
        results.push_back(it->str());
    }

    if (results.empty() || results.size() > 1) { throw std::runtime_error("Error in WURCS Residue Order"); }

    std::string order_string = results[0];
    order_string = order_string.substr(1, order_string.size() - 2);

    std::vector<std::string> order = Utils::split(order_string, '-');
    std::vector<int> residue_order;
    for (auto &it: order) {
        residue_order.emplace_back(atoi(it.c_str()));
    }
    return residue_order;
}

std::vector<std::string> Sails::WURCS::extract_wurcs_linkage_order(const std::string &wurcs) {
    std::regex regex("[a-z]{1}[0-9]{1}-{1}[a-z]{1}[0-9]");
    std::vector<std::string> results;
    auto matches_begin = std::sregex_iterator(wurcs.begin(), wurcs.end(), regex);
    auto matches_end = std::sregex_iterator();

    for (auto it = matches_begin; it != matches_end; ++it) {
        results.push_back(it->str());
    }

    return results;
}


std::vector<std::string> Sails::WURCS::form_residue_name_order(ResidueDatabase &residue_database,
                                                               std::vector<std::string> unique_residues,
                                                               std::vector<int> residue_order) {
    std::vector<std::string> unique_residue_names;
    for (auto &wurcs_residue_id: unique_residues) {
        std::optional<std::string> key = get_key(residue_database, [&](const ResidueData &data) {
            return data.wurcs_code == wurcs_residue_id;
        });
        if (!key.has_value()) {
            throw std::runtime_error("Cannot create psuedo-glycan, an unknown residue was found.");
        }
        unique_residue_names.emplace_back(key.value());
    }

    std::vector<std::string> residue_name_order;
    residue_name_order.reserve(residue_order.size());
    for (auto &order: residue_order) {
        residue_name_order.emplace_back(unique_residue_names[order - 1]);
    }

    return residue_name_order;
}

gemmi::Structure Sails::WURCS::generate_pseudo_structure() {
    gemmi::Structure pseudo_structure;
    gemmi::Model pseudo_model;
    gemmi::Chain chain = gemmi::Chain("A");

    pseudo_model.chains.emplace_back(chain);
    pseudo_structure.models.emplace_back(pseudo_model);
    return pseudo_structure;
}

void Sails::WURCS::add_linkage_to_pseudo_glycan(std::vector<Glycosite> &glycosites, PseudoGlycan &pseudo_glycan,
                                                const std::string &linkage_string) {
    std::vector<std::string> split_linkage = Utils::split(linkage_string, '-');
    const std::string &donor_linkage = split_linkage[0];
    const std::string &acceptor_linkage = split_linkage[1];

    int donor_sugar = donor_linkage[0] - 'a';
    std::string donor_atom = "C";
    donor_atom += donor_linkage[1];
    int acceptor_sugar = acceptor_linkage[0] - 'a';
    std::string acceptor_atom = "C";
    acceptor_atom += acceptor_linkage[1];

    pseudo_glycan.add_linkage(glycosites[donor_sugar], glycosites[acceptor_sugar], donor_atom, acceptor_atom);
}

//WURCS=2.0/3,7,6/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3-3-3/a4-b1_b4-c1_c3-d1_c6-e1_e3-f1_f2-g1
Sails::PseudoGlycan Sails::WURCS::generate_pseudo_glycan(const std::string &wurcs,
                                                   gemmi::Structure *structure,
                                                   Glycosite &glycosite,
                                                   LinkageDatabase &linkage_database,
                                                   ResidueDatabase &residue_database) {
    std::vector<std::string> unique_residues = extract_wurcs_unique_residues(wurcs);
    std::vector<int> residue_order = extract_wurcs_residue_order(wurcs);
    std::vector<std::string> linkage_order = extract_wurcs_linkage_order(wurcs);

    std::vector<std::string> residue_name_order = form_residue_name_order(
        residue_database, unique_residues, residue_order);

    if (structure != nullptr) {
        structure->models[0].chains.emplace_back("X");
    }
        // gemmi::Structure pseudo_structure = generate_pseudo_structure();
        // structure = &pseudo_structure;

    std::vector<Glycosite> glycosites;
    for (int i = 0; i < residue_name_order.size(); i++) {
        gemmi::Residue residue;
        residue.seqid = gemmi::SeqId(std::to_string(i));
        residue.name = residue_name_order[i];
        auto last_chain = Utils::get_last_chain(structure);
        if (!last_chain.has_value()) {throw std::runtime_error("The structure must contain at least one chain.");}
        last_chain.value()->residues.emplace_back(residue);
        int last_chain_index = Utils::get_last_chain_index(structure);
        glycosites.emplace_back(0, last_chain_index, i);
    }

    PseudoGlycan pseudo_glycan = {structure, residue_database, glycosite};
    pseudo_glycan.add_sugar("", 0, glycosite);

    for (int i = 0; i < residue_name_order.size(); i++) {
        pseudo_glycan.add_sugar("", i+1, glycosites[i]);
    }

    pseudo_glycan.add_linkage(glycosite, glycosites[0], "ND2","C1");
    for (const auto &linkage: linkage_order) {
        add_linkage_to_pseudo_glycan(glycosites, pseudo_glycan, linkage);
    }

    return pseudo_glycan;
}

template<typename key, typename value, typename Predicate>
std::optional<key> Sails::WURCS::get_key(const std::map<key, value> &map, Predicate predicate) {
    for (const auto &[k, v]: map) {
        if (predicate(v)) {
            return k;
        }
    }
    return std::nullopt;
}
