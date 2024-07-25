//
// Created by Jordan Dialpuri on 24/07/2024.
//

#ifndef SAILS_DOT_H
#define SAILS_DOT_H

#include <stdio.h>
#include <string>
#include "sails-glycan.h"
#include "sails-topology.h"
#include "sails-json.h"
#include "sails-sequence.h"

namespace Sails {

    struct Dot {
        explicit Dot(gemmi::Structure& structure);

        std::map<Glycosite, std::string> get_all_dotfiles();

        std::string get_dotfile(Sails::Glycosite& glycosite);

    private:
        static std::string header();
        static std::string footer();
        std::string get_format(const Glycosite &site);

    private:
        gemmi::Structure m_structure;
        ResidueDatabase m_database;
        Topology m_topology;
    };
}

#endif //SAILS_DOT_H
