//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_UTILS_H
#define SAILS_SAILS_UTILS_H

#include <iostream>

#include "gemmi/math.hpp"
#include "gemmi/model.hpp"

#include "sails-model.h"

namespace Sails::Utils {
    inline void print_vector(const gemmi::Vec3& vec) {
        std::cout << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
    }

    inline std::string format_residue_key(const gemmi::Residue* residue) {
        return residue->name + "-" + residue->seqid.str();
    }

    gemmi::Chain get_chain_from_glycosite(const Glycosite& site, const gemmi::Structure& structure) {
        return structure.models[site.model_idx].chains[site.chain_idx];
    }

    gemmi::Residue get_residue_from_glycosite(const Glycosite& site,const gemmi::Structure& structure) {
        return structure.models[site.model_idx].chains[site.chain_idx].residues[site.residue_idx];
    }

    gemmi::Atom get_atom_from_glycosite(const Glycosite& site, const gemmi::Structure& structure) {
        if (site.atom_idx == -1) {throw std::runtime_error("Site has not been initialised from a Mark");}
        return structure.models[site.model_idx].chains[site.chain_idx].residues[site.residue_idx].atoms[site.atom_idx];
    }

} // namespace Sails::Utils


#endif //SAILS_SAILS_UTILS_H
