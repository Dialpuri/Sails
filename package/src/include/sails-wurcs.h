//
// Created by jordan on 05/11/24.
//

#ifndef SAILS_SAILS_WURCS_H
#define SAILS_SAILS_WURCS_H

#include "sails-glycan.h"

namespace Sails {

    class WURCS {
    public:
        WURCS() = default;

        static std::string generate_wurcs(Sails::Glycan* glycan, Sails::ResidueDatabase& residue_database);

        static Sails::Glycan generate_psuedo_glycan(const std::string& wurcs, gemmi::Structure* structure);

    private:
        static std::string get_unit_count(Sails::Glycan* glycan);

        static std::string get_unique_residue_list(Sails::Glycan* glycan, Sails::ResidueDatabase& residue_database);

        static std::string get_residue_order(Sails::Glycan* glycan, Sails::ResidueDatabase& residue_database);

        static std::string get_wurcs_header() {return "WURCS=2.0";}

        static std::string get_section_delimiter() {return "/";}
    };

}

#endif //SAILS_SAILS_WURCS_H
