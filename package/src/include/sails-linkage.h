//
// Created by Jordan Dialpuri on 07/07/2024.
//

#ifndef SAILS_SAILS_LINKAGE_H
#define SAILS_SAILS_LINKAGE_H

#include "sails-glycan.h"
#include "sails-utils.h"
#include "sails-vector.h"
#include <gemmi/cif.hpp>             // for cif::read_file
#include <string>


namespace Sails {
    struct Model {
        Model() = default;
        Model(gemmi::Structure &structure, LinkageDatabase &linkage_database, ResidueDatabase &residue_database)
                : structure(structure), linkage_database(linkage_database), residue_database(residue_database) {
            monomer_library_path = Sails::Utils::get_environment_variable("CLIBD") + "/monomers";
        }

        Sails::Glycan extend(Sails::Glycan &glycan, int base_seqid);

    private:

        std::optional<gemmi::Residue> get_monomer(const std::string& monomer);

        gemmi::Residue calculate_residue_translation(const gemmi::Residue& residue, LinkageData& data);

        std::vector<gemmi::Atom *> extract_atoms(const gemmi::Residue &residue, LinkageData &data, std::optional<gemmi::Residue> library_monomer, std::vector<gemmi::Atom *> reference_atoms);

    private:
        gemmi::Structure structure;
        Sails::LinkageDatabase linkage_database;
        Sails::ResidueDatabase residue_database;
        std::string monomer_library_path = "";
    };


}

#endif //SAILS_SAILS_LINKAGE_H
