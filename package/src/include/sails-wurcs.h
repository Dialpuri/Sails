//
// Created by jordan on 05/11/24.
//

#ifndef SAILS_SAILS_WURCS_H
#define SAILS_SAILS_WURCS_H

#include "sails-glycan.h"
#include <regex>

namespace Sails {
    class WURCS {
    public:
        WURCS() = default;

        static std::string generate_wurcs(Glycan *glycan, ResidueDatabase &residue_database);

        static PseudoGlycan generate_pseudo_glycan(const std::string &wurcs, gemmi::Structure *structure,
                                             Glycosite &glycosite,
                                             LinkageDatabase &linkage_database, ResidueDatabase &residue_database);

    private:
        static std::string get_unit_count(Sails::Glycan *glycan);

        static std::string get_unique_residue_list(Sails::Glycan *glycan, Sails::ResidueDatabase &residue_database);

        static std::string get_residue_order(Sails::Glycan *glycan, Sails::ResidueDatabase &residue_database);

        static std::string get_link_list(Sails::Glycan *glycan, Sails::ResidueDatabase &residue_database);

        static std::string get_wurcs_header() { return "WURCS=2.0"; }

        static std::string get_section_delimiter() { return "/"; }

        static std::vector<int> calculate_residue_order(Glycan *glycan, ResidueDatabase &residue_database);

        // PARSING

        static std::vector<std::string> extract_wurcs_unique_residues(const std::string &wurcs);

        static std::vector<int> extract_wurcs_residue_order(const std::string &wurcs);

        static std::vector<std::string> extract_wurcs_linkage_order(const std::string &wurcs);

        static std::vector<std::string> form_residue_name_order(ResidueDatabase &residue_database,
                                                               std::vector<std::string> unique_residues,
                                                               std::vector<int> residue_order);

        static gemmi::Structure generate_pseudo_structure(
        );

        static void add_linkage_to_pseudo_glycan(std::vector<Sails::Glycosite> &glycosites,
                                                 Sails::PseudoGlycan &pseudo_glycan,
                                                 const std::string &linkage_string);

        /**
         * @brief Get key from map using predicate, used to query a linkage or residue database
         *
         * @param map The map to get a key from
         * @param predicate Predicate to search using
         *
         * @return optional key
         */
        template<typename key, typename value, typename Predicate>
        static std::optional<key> get_key(const std::map<key, value> &map, Predicate predicate);
    };
}

#endif //SAILS_SAILS_WURCS_H
