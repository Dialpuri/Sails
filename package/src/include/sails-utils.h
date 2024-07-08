//
// Created by jordan on 05/07/24.
//
#pragma once
#ifndef SAILS_SAILS_UTILS_H
#define SAILS_SAILS_UTILS_H

#include <iostream>
#include <fstream>

#include "gemmi/math.hpp"
#include "gemmi/model.hpp"

#include "sails-model.h"

namespace Sails::Utils {
    std::string get_environment_variable(std::string const &key);

    void print_vector(const gemmi::Vec3 &vec);

    void print_position(const gemmi::Position &position);

    std::string format_residue_key(const gemmi::Residue *residue);

    gemmi::Chain get_chain_from_glycosite(const Glycosite &site, const gemmi::Structure &structure);

    gemmi::Residue get_residue_from_glycosite(const Glycosite &site, const gemmi::Structure &structure);

    gemmi::Atom get_atom_from_glycosite(const Glycosite &site, const gemmi::Structure &structure);

    std::string linkage_to_id(const Sails::LinkageData &data);

    void save_residues_to_file(std::vector<gemmi::Residue>& residues, const std::string& path);

} // namespace Sails::Utils


#endif //SAILS_SAILS_UTILS_H
