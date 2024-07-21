//
// Created by Jordan Dialpuri on 06/07/2024.
//

#ifndef SAILS_SAILS_SEQUENCE_H
#define SAILS_SAILS_SEQUENCE_H

#include <iostream>

#include "gemmi/model.hpp" // for structure
#include "gemmi/resinfo.hpp" // for tabulated_residue

#include "sails-model.h"

namespace Sails {
    /**
     * @brief Finds N-glycosylation sites in a given structure.
     *
     * This function searches for N-glycosylation sites in a given structure. An N-glycosylation site is defined as a
     * sequence of three consecutive residues: N-X-S/T, where N is the first residue, S/T is the third residue, and X
     * can be any residue except proline (P). Only sites that satisfy this sequence pattern are considered as
     * N-glycosylation sites. The function returns a vector of Glycosite objects that represent the position of each
     * N-glycosylation site in the structure.
     *
     * @param structure The structure in which to search for N-glycosylation sites.
     * @return A vector of Glycosite objects representing the position of each N-glycosylation site in the structure.
     *         If no N-glycosylation sites are found, an empty vector is returned.
     */
    Glycosites find_n_glycosylation_sites(const gemmi::Structure &structure);

    /**
     * @brief Finds C-glycosylation sites in a given structure.
     *
     * This function searches for C-glycosylation sites in a given structure. An C-glycosylation site is defined as a
     * sequence of three consecutive residues: W-X-X-W, where either W can be mannosylated.
     * Only sites that satisfy this sequence pattern are considered as C-glycosylation sites. The function returns a
     * vector of Glycosite objects that represent the position of each C-glycosylation site in the structure.
     *
     * @param structure The structure in which to search for C-glycosylation sites.
     * @return A vector of Glycosite objects representing the position of each C-glycosylation site in the structure.
     *         If no C-glycosylation sites are found, an empty vector is returned.
     */
    Glycosites find_c_glycosylation_sites(const gemmi::Structure &structure);
}

#endif //SAILS_SAILS_SEQUENCE_H
