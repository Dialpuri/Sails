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
#include "sails-cif.h"

#include <string>

#include <gemmi/cif.hpp>             // for cif::read_file
#include <gemmi/modify.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/to_pdb.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_cif.hpp>


namespace Sails {

 /**
   * @struct SuperpositionResult
   * @brief Represents the result of a superposition operation between two residues.
   *
   * @param new_residue The newly transformed residue.
   * @param transformation The transformation matrix used for superposition.
   * @param reference_residue The reference residue used for superposition.
   *
   */
 struct SuperpositionResult {

  SuperpositionResult() = default;

  SuperpositionResult(gemmi::Residue new_residue, const gemmi::Transform &transformation,
                      gemmi::Residue reference_residue)
   : new_residue(std::move(new_residue)),
     transformation(transformation),
     reference_residue(std::move(reference_residue)) {
  }

  /**
   * @brief Checks if the transformation vector is null.
   *
   * @return True if the transformation vector has NaN values, otherwise False.
   */
  [[nodiscard]] bool is_null() const { return transformation.vec.has_nan(); }

  gemmi::Residue new_residue;
  gemmi::Transform transformation;
  gemmi::Residue reference_residue;
 };

 /**
  * @class Model
  * @brief Represents a model used for extending glycans and performing other operations.
  */
 struct Model {
  friend struct TorsionAngleRefiner;
  Model() = default;

  Model(gemmi::Structure *structure, LinkageDatabase &linkage_database, ResidueDatabase &residue_database)
   : structure(structure), linkage_database(linkage_database), residue_database(residue_database) {
   monomer_library_path = Utils::get_environment_variable("CLIBD") + "/monomers";
  }

  /**
   * Extends the given glycan by adding new sugars based on the linkage database.
   *
   * @param glycan The glycan to be extended.
   * @param base_glycosite The base sequence ID to determine terminal sugars.
   * @param density The density class used to score residues into experimental data.
   * @return The extended glycan.
   */
  Glycan extend(Glycan &glycan, Glycosite &base_glycosite, Density &density, bool debug);


  /**
   * @brief Save the model with link records to a specified path.
   *
   * This method saves the model with links links to a specified path. It creates a CIF document
   * using the given structure and adds the links to the "_struct_conn" category.
   * The CIF document is then written to the specified path.
   *
   * @param path The path where the CIF file should be saved.
   * @param links A vector of LinkRecord objects representing the links to be saved.
   *              Each LinkRecord contains the labels for a link.
   */
  void save(const std::string &path, std::vector<LinkRecord> &links);

  /**
   * @brief Saves the model to a file.
   *
   * This method saves the model to the specified file path. The method first opens
   * a new output file stream using the given path. Then, it writes the model data in
   * CIF format to the output stream using the `write_cif_to_stream` function. Finally,
   * it closes the output file stream.
   *
   * @param path The file path where the model will be saved.
   */
  void save(const std::string &path);

  /**
   * @brief Get the structure.
   *
   * This method returns a pointer to the structure object.
   *
   * @return A pointer to the structure object.
   */
  [[nodiscard]] gemmi::Structure *get_structure() const { return structure; }

 private:
  typedef std::map<int, std::vector<Sails::SuperpositionResult> > PossibleAdditions;

  /**
   * Extends the given glycan by adding new sugars based on the linkage database.
   *
   * @param density The density class used to score residues into experimental data.
   * @param debug A flag indicating whether debug output should be enabled.
   * @param is_sugar_only_chain A boolean reference that will be updated to indicate whether the extended glycan is a sugar-only chain.
   * @param terminal_sugar The terminal sugar from which to extend the glycan.
   */
  void extend_if_possible(Density &density, bool debug, bool &is_sugar_only_chain,
                          const Sugar *terminal_sugar);

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

  /**
   * Finds favoured additions based on a terminal sugar and a list of possible additions.
   * A favoured addition is one where the depth of the new residue matches one of the preferred depths in the
   * residue database for that residue name.
   *
   * @param terminal_sugar The terminal sugar to consider.
   * @param possible_additions The list of possible additions for the given terminal sugar.
   * @return A vector of SuperpositionResult objects representing the favoured additions.
   */
  std::vector<SuperpositionResult> find_favoured_additions(const Sugar *&terminal_sugar, PossibleAdditions &
                                                           possible_additions);

  /**
   * @brief Adds a sugar to a structure.
   *
   * This method adds a sugar to a structure given the terminal sugar, superposition result,
   * and a flag indicating whether the sugar chain is only composed of sugars.
   *
   * @param terminal_sugar The terminal sugar to add.
   * @param result The superposition result representing the sugar to be added.
   * @param is_sugar_only_chain A flag indicating whether the sugar chain is only composed of sugars.
   */
  void add_sugar_to_structure(const Sugar *terminal_sugar, SuperpositionResult &result, bool is_sugar_only_chain);

  /**
   * @brief Calculates the clash score for the given SuperpositionResult.
   *
   * The clash score is calculated by finding the number of nearby atoms for each atom in the
   * SuperpositionResult. Nearby atoms are found using a NeighborSearch with a given radius.
   *
   * @param result The SuperpositionResult from which to calculate the clash score.
   * @return The calculated clash score.
   */
  [[nodiscard]] double calculate_clash_score(const SuperpositionResult &result) const;


  /**
   * @brief Removes the leaving atom from the given residue objects.
   *
   * The leaving atom is identified by the acceptor number in the linkage data.
   * The leaving atom is represented as "O" followed by the acceptor number.
   * The leaving atom is removed from both the new_monomer and reference_library_monomer residue objects.
   *
   * @param data The linkage data that contains the acceptor number.
   * @param reference_library_monomer The reference library monomer residue object.
   * @param new_monomer The new monomer residue object.
   */
  static void remove_leaving_atom(LinkageData &data, gemmi::Residue &reference_library_monomer,
                                  gemmi::Residue &new_monomer);

  /**
   * @brief Move the positions of acceptor atoms based on given parameters.
   *
   * This method calculates a new position for each acceptor atom in the provided vector
   * of atoms, based on the length, angles, and torsions given. The new position is then
   * assigned to the next atom in the vector.
   *
   * @param atoms The vector of atoms representing the acceptor atoms.
   * @param length The length parameter used for calculating the new positions.
   * @param angles The vector of angles used for calculating the new positions.
   * @param torsions The vector of torsions used for calculating the new positions.
   */
  static void move_acceptor_atomic_positions(std::vector<gemmi::Atom *> &atoms, double length,
                                             std::vector<double> &angles, std::vector<double> &torsions);

  /**
   * @brief Applies a superposition transformation to a set of atoms.
   *
   * This method calculates the superposition transformation between a set of atoms and a set of reference atoms,
   * specified by their positions. The transformation consists of a translation and rotation applied to the set of atoms,
   * such that they align as closely as possible with the reference atoms.
   *
   * @param atoms The set of atoms to be transformed.
   * @param reference_atoms The set of reference atoms used for the superposition calculation.
   * @param length The length parameter used in the superposition calculation.
   * @param angles The angles parameter used in the superposition calculation.
   * @param torsions The torsions parameter used in the superposition calculation.
   * @return The calculated transformation representing the superposition.
   */
  static gemmi::Transform superpose_atoms(std::vector<gemmi::Atom *> &atoms,
                                          std::vector<gemmi::Atom> &reference_atoms, double length,
                                          std::vector<double> &angles, std::vector<
                                           double> &torsions);

  /**
     * Retrieves the monomer with the given name from the monomer library.
     *
     * @param monomer The name of the monomer to retrieve.
     * @param remove_h
     * @return An optional value that contains the monomer, or an empty optional if the monomer does not exist.
     */
  std::optional<gemmi::Residue> get_monomer(const std::string &monomer, bool remove_h);


  /**
   * @brief Rotates the exocyclic atoms of a residue and finds the best position based on density scoring.
   *
   * This method performs a rotation of the exocyclic atoms of a given residue around an axis defined by two atoms.
   * It calculates the rotation matrix and applies it to the exocyclic atom specified by the third atom.
   * The resulting positions are evaluated based on density scoring, and the best position is assigned to the exocyclic atom.
   *
   * @param residue Pointer to the gemmi::Residue object representing the residue.
   * @param atoms Vector of strings representing the names of the three atoms that define the rotation axis.
   * @param density Reference to a Density object used for scoring the rotated positions.
   */
  static void rotate_exocyclic_atoms(gemmi::Residue *residue, std::vector<std::string> &atoms, Density &density);

  /**
   * @brief Checks if a chain consists only of sugars
   *
   */
  bool check_if_sugar_only_chain(std::vector<Sugar *> sugars);

  /**
  * @brief Checks if a residue is in the residue database
  *
  */
  bool check_entry_in_database(gemmi::Residue *residue);


  /**
   * @brief Print the addition log of a new residue to a terminal sugar.
   *
   * This function prints the addition log by outputting the new residue being added to the terminal sugar.
   *
   * @param terminal_sugar The terminal sugar where the new residue is being added.
   * @param addition The superposition result containing the new residue being added.
   */
  void print_addition_log(const Sugar *terminal_sugar, SuperpositionResult &addition);

  /**
   * @brief Print the log for an attempted addition of a residue to a model.
   *
   * This function prints a log message indicating the attempted addition of a residue to a model.
   * The log message includes information about the donor and acceptor residues, as well as the donor and acceptor numbers.
   *
   * @param residue_ptr A pointer to the gemmi::Residue object representing the residue being added.
   * @param data The LinkageData object containing information about the donor and acceptor residues.
   */
  void print_attempted_addition_log(gemmi::Residue *residue_ptr, LinkageData &data);

  /**
   * @brief Prints the rejection log.
   *
   * This method prints the rejection log to the standard output.
   */
  void print_rejection_log();

  /**
   * @brief Prints a successful log message after adding a residue to the density.
   *
   * This method takes in a Density object and an optional SuperpositionResult object as parameters.
   * It calculates the RSCC score using the SuperpositionResult and prints a log message indicating the
   * successful addition of the residue along with the RSCC score.
   *
   * @param density The Density object representing the density in which the residue is added.
   * @param opt_result (optional) The SuperpositionResult object containing the new residue and other details.
   *
   * @note If the opt_result parameter is not provided, the method will throw an exception.
   */
  void print_successful_log(Density &density, std::optional<SuperpositionResult> opt_result);

 private:
  gemmi::Structure *structure{};
  LinkageDatabase linkage_database;
  ResidueDatabase residue_database;
  std::string monomer_library_path;
  std::string special_monomer_path = "package/models/special_monomers";
  std::map<std::string, gemmi::Residue> monomers;
 };
}

#endif //SAILS_SAILS_LINKAGE_H
