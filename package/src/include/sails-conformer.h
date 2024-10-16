//
// Created by Jordan Dialpuri on 22/09/2024.
//

#ifndef SAILS_CONFORMER_H
#define SAILS_CONFORMER_H

#include "sails-linkage.h"
#include "sails-model.h"

namespace Sails {
    struct Conformer : private Model {
        Conformer() = default;

        Conformer(ResidueDatabase &residue_database, LinkageDatabase &linkage_database): Model(
            nullptr, linkage_database, residue_database) {
        }

        void superimpose(const std::string &residue_1_name, const std::string &residue_2_name,
                         int donor_number, int acceptor_number, double phi, double psi);


        void get_all_conformers(const std::string &residue_1_name, const std::string &residue_2_name, int donor_number,
                                int acceptor_number, int step);

        void remove_leaving_atoms(Sails::LinkageData &data, gemmi::Residue &reference_library_monomer,
                                       gemmi::Residue &new_monomer);
    };
}

#endif //SAILS_CONFORMER_H
