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
            monomer_library_path = Utils::get_environment_variable("CLIBD") + "/monomers";
        }

        /**
         * Extends the given glycan by adding new sugars based on the linkage database.
         *
         * @param glycan The glycan to be extended.
         * @param base_seqid The base sequence ID to determine terminal sugars.
         * @return The extended glycan.
         */
        Glycan extend(Glycan &glycan, int base_seqid);

    private:
        /**
         * Retrieves the monomer with the given name from the monomer library.
         *
         * @param monomer The name of the monomer to retrieve.
         * @return An optional value that contains the monomer, or an empty optional if the monomer does not exist.
         */
        std::optional<gemmi::Residue> get_monomer(const std::string& monomer);

        /**
         * Performs the translation of a residue based on the given linkage data.
         *
         * The donor atoms of the input residue are used to calculate a superposition of a given new monomer. A new
         * monomer is loaded and transformed to the correct position and returned.
         *
         * @param residue The residue to calculate the translation for.
         * @param data The linkage data containing donor and acceptor information.
         * @return The translated residue.
         * @throws std::runtime_error if any required atom or data is not found, or unexpected atom count.
         */
        gemmi::Residue calculate_residue_translation(const gemmi::Residue& residue, LinkageData& data);

    private:
        gemmi::Structure structure;
        LinkageDatabase linkage_database;
        ResidueDatabase residue_database;
        std::string monomer_library_path;
    };


}

#endif //SAILS_SAILS_LINKAGE_H
