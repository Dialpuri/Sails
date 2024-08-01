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
 /**
  * @struct Dot
  *
  * @brief The Dot class is used to generate dotfiles for glycosites in a given structure.
  *
  * This class generates dotfiles for glycosites in a given structure, which can then be used to create an SNFG.
  * The dotfiles can be obtained for all glycosites or for a specific glycosite.
  *
  * The class also provides methods to add a header and footer to the dotfile, as well as
  * a method to identify the format of a glycosite.
  *
  */
 struct Dot {
  /**
   * @brief Constructor for Dot object.
   *
   * @param structure The structure object to use for generating Dot files, which are used to create SNFGs.
   *
   * @return None.
   */
  explicit Dot(gemmi::Structure& structure);

  /**
   * @brief Retrieves all dot files for all N-glycosylation glycosites in the structure.
   *
   * This method returns a map where the Glycosite is the key and the associated
   * dot file is the value. The dot file contains information about the glycosite
   * and its surrounding residues in a specific structure.
   *
   * @returns A map where the Glycosite is the key and the associated dot file is the value.
   */
  std::map<Glycosite, std::string> get_all_dotfiles();

  /**
   * @brief Generates a DOT file representation of the given glycosite.
   *
   * This method takes a Glycosite object and returns a string that represents a DOT file
   * containing the topology information of the glycosite.
   *
   * @param glycosite The glycosite for which the DOT file will be generated.
   * @return A string representing the DOT file of the glycosite.
   */
  std::string get_dotfile(Sails::Glycosite& glycosite);

 private:

  /**
   * @brief Return the header string for the Glycan graph.
   *
   * This method returns the header string for the Glycan graph in DOT format.
   * The header string includes properties such as rank direction, padding, node separation,
   * and rankset.
   *
   * @return The header string for the Glycan graph.
   */
  static std::string header();

  /**
   * @brief Generates a footer for the dot file.
   *
   * This method returns a string that represents the footer of a dot file.
   *
   * @return The footer string.
   */
  static std::string footer();

  /**
   * @brief Get the format string for the given glycosite.
   *
   * This method retrieves the format string for the given glycosite.
   * The format string contains information about the fill color,
   * shape, and label of the glycosite.
   *
   * @param site The glycosite for which the format is requested.
   *
   * @return The format string for the given glycosite.
   */
  std::string get_format(const Glycosite &site);

  /**
   * Internal structure - used to find linkages
   */
 private:
        gemmi::Structure m_structure;
  /**
   * Internal residue database - used to find linkages
   */
  ResidueDatabase m_database;
  /**
   * Internal topology - used to find linkages
   */
  Topology m_topology;
    };
}

#endif //SAILS_DOT_H
