//
// Created by Jordan Dialpuri on 08/08/2024.
//

#ifndef SAILS_SOLVENT_H
#define SAILS_SOLVENT_H

#include <iostream>

#include <gemmi/grid.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mtz.hpp>
#include <gemmi/solmask.hpp>

#include "sails-model.h"
#include "sails-utils.h"

namespace Sails {
 /**
  * @class SolventAccessibility
  * @brief Class for calculating solvent accessibility of glycosites in a structure.
  *
  * This class provides methods to calculate the solvent accessibility of glycosites
  * in a given structure. The solvent accessibility is calculated using a grid-based
  * approach where if a given residue is within n Angstroms of the surface of the solvent mask, it
  * is solvent accessible
  */
 class SolventAccessibility {
 public:
  typedef std::map<Glycosite, double> SolventAccessibilityMap;


  /**
   * @brief Constructor for SolventAccessibility class.
   *
   * @param structure A pointer to a gemmi::Structure object.
   *
   * This constructor initializes a SolventAccessibility object with the given gemmi::Structure object.
   * It sets up the solvent mask from the structure and calculates the solvent mask.
   */
  explicit SolventAccessibility(gemmi::Structure *structure): m_structure(structure) {
   m_solvent_mask.setup_from(*structure, 1);
   calculate_solvent_mask();
  }

  /**
   * @brief Average and flatten the solvent accessibility map
   *
   * @param map The original solvent accessibility map.
   * @return The flattened solvent accessibility map, where the average solvent accessibility
   *         for each residue site is calculated.
   */
  SolventAccessibilityMap average_solvent_accessibilty_map(const SolventAccessibilityMap &map);

  /**
   * Calculates the solvent accessibility for each glycosite in the structure.
   *
   * @return A map containing the solvent accessibility value for each glycosite.
   */
  std::map<Glycosite, double> calculate_solvent_accessibility();

  /**
   * @brief  Calculates the bounding box of the structure.
   *
   * This method calculates the bounding box of the structure by iterating over
   * all models, chains, residues, and atoms in the structure and extending the
   * box to include each atom position. The calculated box is then expanded by
   * adding a margin of 8 Angstroms to each side.
   *
   * @return  The calculated bounding box of the structure.
   */
 private:
  gemmi::Box<gemmi::Position> calculate_box();

  /**
   * @brief Create a solvent accessibility map.
   *
   * This method creates a solvent accessibility map for the given structure.
   * The solvent accessibility map is a mapping of the glycosite (model index, chain index, residue index, atom index)
   * to its corresponding solvent accessibility value.
   *
   * The solvent accessibility value is initially set to 0 for each glycosite.
   *
   * @return The created solvent accessibilty map.
   */
  SolventAccessibilityMap create_solvent_accessibility_map();

  /**
   * @brief Calculates the solvent mask for the structure.
   *
   * It applies a mask on the grid based on the atomic radii set and the model index of the structure.
   */
  void calculate_solvent_mask();

 private:
  /**
   * @brief Pointer to a gemmi::Structure for use when calculating solvent mask and bounding box
   */
  gemmi::Structure *m_structure;

  /**
   * @brief Solvent mask calculated during constructor
   */
  gemmi::Grid<> m_solvent_mask;
 };
}
#endif //SAILS_SOLVENT_H
