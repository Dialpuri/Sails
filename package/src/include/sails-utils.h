//
// Created by jordan on 05/07/24.
//
#pragma once
#ifndef SAILS_SAILS_UTILS_H
#define SAILS_SAILS_UTILS_H

#include <iostream>
#include <fstream>
#include <gemmi/math.hpp>
#include <gemmi/model.hpp>
#include <gemmi/ccp4.hpp>
#include <gemmi/to_pdb.hpp>

#include "sails-model.h"

namespace Sails::Utils {
 /**
  * @brief Retrieves the value of the specified environment variable.
  *
  * This function retrieves the value of the specified environment variable identified by the provided key.
  *
  * @param key The name of the environment variable to retrieve.
  *
  * @return The value of the specified environment variable. If the environment variable is not found,
  *         an empty string is returned.
  */
 std::string get_environment_variable(std::string const &key);

 /**
  * @brief Prints the components of the given vector.
  *
  * This method prints the x, y, and z components of the provided vector using the format "(x, y, z)" using the
  * std::cout stream
  *
  * @param vec The gemmi::Vec3 object whose coordinates will be printed.
  */
 void print_vector(const gemmi::Vec3 &vec);

 /**
    * @brief Prints the components of the given vector with given label prefix.
    *
    * This method prints the x, y, and z components of the provided vector using the format "(x, y, z)" using the
    * std::cout stream
    *
    * @param vec The gemmi::Vec3 object whose coordinates will be printed.
    */
 void print_vector(const gemmi::Vec3 &vec, const std::string &name);

 /**
  * @brief Prints the coordinates of a given gemmi::Position object.
  *
  * This method prints the x, y, and z components of the provided position using the format "(x, y, z)" using the
  * std::cout stream
  *
  * @param position The gemmi::Position object whose coordinates will be printed.
  */
 void print_position(const gemmi::Position &position);

 /**
  * @brief Formats the key for a given residue.
  *
  * This function formats the key for the provided residue by concatenating its name and sequence identifier.
  *
  * @param residue A pointer to the gemmi::Residue object representing the residue.
  *
  * @return The formatted key for the residue.
  */
 std::string format_residue_key(const gemmi::Residue *residue);


 /**
  * @brief Formats the key for a given glycosite.
  *
  * This function formats the key for the provided residue by concatenating its the model, chain, residue and atom ids.
  *
  * @param glycosite A glycosite
  *
  * @return The formatted key for the site.
  */
 std::string format_site_key(const Glycosite& glycosite);


 /**
  * @brief Retrieves the gemmi::Chain object corresponding to a given Glycosite.
  *
  * This function retrieves the gemmi::Chain object that represents the chain of the
  * provided Glycosite within the specified gemmi::Structure.
  *
  * @param site The glycosite object containing the model, chain, and residue indices.
  * @param structure The gemmi::Structure object from which the chain is to be retrieved.
  *
  * @return The gemmi::Chain object from the specified glycosute.
  */
 gemmi::Chain get_chain_from_glycosite(const Glycosite &site, const gemmi::Structure &structure);

 /**
  * @brief Retrieves the gemmi::Residue corresponding to a given Glycosite.
  *
  * This function retrieves the gemmi::Residue object that represents the residue of the
  * provided Glycosite within the specified gemmi::Structure.
  *
  * @param site The glycosite object containing the model, chain, and residue indices.
  * @param structure The structure object from which to retrieve the residue.
  *
  * @return The gemmi::Residue from the specified glycosite.
  */
 gemmi::Residue get_residue_from_glycosite(const Glycosite &site, const gemmi::Structure &structure);

 /**
  * @brief Retrieves the gemmi::Atom from the specified glycosite.
  *
  * This function retrieves the gemmi::Atom object that represents the atom of the
  * provided Glycosite within the specified gemmi::Structure.
  *
  * @param site The glycosite object containing the model, chain, and residue indices.
  * @param structure The structure from which to retrieve the atom.
  *
  * @return The gemmi::Atom from the specified glycosite.
  * @throws std::runtime_error if the site has not been initialised from a Mark.
  */
 gemmi::Atom get_atom_from_glycosite(const Glycosite &site, const gemmi::Structure &structure);

 /**
  * @brief Converts the given LinkageData object into an ID string.
  *
  * This function takes a reference to a LinkageData object and converts it into an ID string. The ID string is constructed
  * by concatenating the donor, donor_number, acceptor_number, and acceptor fields of the LinkageData object.
  *
  * @param data The LinkageData object to convert to an ID string.
  *
  * @return The ID string representation of the LinkageData object.
  */
 std::string linkage_to_id(const Sails::LinkageData &data);

 /**
  * @brief Save a vector of gemmi::Residue objects to a specified file path.
  *
  * This method takes a vector of gemmi::Residue objects and saves them to the specified file path.
  * It creates a gemmi::Structure, gemmi::Model, and gemmi::Chain to hold the residues.
  * The provided residues are appended to the chain, and the chain is added to the model.
  * Finally, the model is added to the structure.
  * The structure is then written to the specified file path using the gemmi::write_pdb function.
  *
  * @param residues The vector of gemmi::Residue objects to save.
  * @param path The file path to save the residues to.
  *
  * @note The file path must be valid and accessible for writing. If the file already exists, it will be overwritten.
  */
 void save_residues_to_file(std::vector<gemmi::Residue> residues, const std::string &path);


 /**
  * @brief Saves a grid to a file in CCP4 format.
  *
  * This function saves the given grid to a file in CCP4 format. The grid is defined by the provided
  * gemmi::Grid<> object and the file path is specified using the std::string path parameter.
  *
  * @param grid The grid object to save to file.
  * @param path The path to the file where the grid will be saved.
  */
 void save_grid_to_file(const gemmi::Grid<> &grid, const std::string &path);
} // namespace Sails::Utils


#endif //SAILS_SAILS_UTILS_H
