//
// Created by Jordan Dialpuri on 07/07/2024.
//

#ifndef SAILS_SAILS_LINKAGE_H
#define SAILS_SAILS_LINKAGE_H

#include "sails-glycan.h"
#include "sails-utils.h"
#include "sails-topology.h"
#include "sails-vector.h"
#include "sails-density.h"

#include <string>

#include <gemmi/cif.hpp>             // for cif::read_file
#include <gemmi/modify.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/to_pdb.hpp>

namespace Sails {


    struct SuperpositionResult {
        SuperpositionResult() = default;
        SuperpositionResult(gemmi::Residue new_residue, const gemmi::Transform &transformation, gemmi::Residue reference_residue)
            : new_residue(std::move(new_residue)),
              transformation(transformation),
              reference_residue(std::move(reference_residue)) {
        }

        bool is_null() const {return transformation.vec.has_nan();}

        gemmi::Residue new_residue;
        gemmi::Transform transformation;
        gemmi::Residue reference_residue;
    };

    struct Model {
        Model() = default;
        Model(gemmi::Structure &structure, LinkageDatabase &linkage_database, ResidueDatabase &residue_database)
                : structure(structure), linkage_database(linkage_database), residue_database(residue_database) {
            monomer_library_path = Utils::get_environment_variable("CLIBD") + "/monomers";
        }

        bool check_entry_in_database(gemmi::Residue residue);

        /**
         * Extends the given glycan by adding new sugars based on the linkage database.
         *
         * @param glycan The glycan to be extended.
         * @param base_seqid The base sequence ID to determine terminal sugars.
         * @param density The density class used to score residues into experimental data.
         * @return The extended glycan.
         */
        Glycan extend(Glycan &glycan, int base_seqid, Density& density);

        static void move_acceptor_atomic_positions(std::vector<gemmi::Atom *>& atoms, double length,
                                                   std::vector<double>& angles, std::vector<double>& torsions);

        static gemmi::Transform superpose_atoms(std::vector<gemmi::Atom *>& atoms, std::vector<gemmi::Atom>& reference_atoms, double length, std::vector<double> &angles, std::vector<
                                                    double> &torsions);

        void save_pdb(const std::string& path);

        void remove_leaving_atom(Sails::LinkageData &data, gemmi::Residue& reference_library_monomer,
                                 gemmi::Residue& new_monomer);

        [[nodiscard]] gemmi::Structure get_structure() const {return structure;}

    private:
        /**
         * Retrieves the monomer with the given name from the monomer library.
         *
         * @param monomer The name of the monomer to retrieve.
         * @param remove_h
         * @return An optional value that contains the monomer, or an empty optional if the monomer does not exist.
         */
        std::optional<gemmi::Residue> get_monomer(const std::string& monomer, bool remove_h);

        /**
         * Performs the translation of a residue based on the given linkage data.
         *
         * The donor atoms of the input residue are used to calculate a superposition of a given new monomer. A new
         * monomer is loaded and transformed to the correct position and returned.
         *
         * @param residue The residue to calculate the translation for.
         * @param data The linkage data containing donor and acceptor information.
         * @param density The constructed density object for accessing experimental data.
         * @param refine Boolean to run real space simplex refinement on the residue.
         * @return Optionally, the translated residue.
         * @throws std::runtime_error if any required atom or data is not found, or unexpected atom count.
         */
        std::optional<Sails::SuperpositionResult> add_residue(
            gemmi::Residue *residue, LinkageData &data, Density &density, bool refine);

        void add_sugar_to_structure(const Sugar* terminal_sugar, SuperpositionResult& result);

        static void rotate_exocyclic_atoms(gemmi::Residue *residue, std::vector<std::string>& atoms, Density &density);

    private:
        gemmi::Structure structure;
        LinkageDatabase linkage_database;
        ResidueDatabase residue_database;
        std::string monomer_library_path;
        std::map<std::string, gemmi::Residue> monomers;
    };


}

#endif //SAILS_SAILS_LINKAGE_H
