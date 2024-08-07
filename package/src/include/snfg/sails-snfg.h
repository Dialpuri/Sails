//
// Created by Jordan Dialpuri on 04/08/2024.
//

#ifndef SAILS_SNFG_H
#define SAILS_SNFG_H

#include <gemmi/model.hpp>

#include "../sails-model.h"
#include "../sails-topology.h"

namespace Sails {
 /**
  * @brief Represents different types of SVG elements.
  *
  * The SVGType enumeration represents different types of SVG (Scalable Vector Graphics) elements that can be used
  * to create visual elements in an SVG document.
  *
  * This enumeration defines the following types:
  * - square: Represents a square element.
  * - circle: Represents a circle element.
  * - line: Represents a line element.
  * - text: Represents a text element.
  */
 enum SVGType {
  square, circle, line, text
 };

 /**
  * @brief Represents an SVG (Scalable Vector Graphics) object.
  *
  * The SVGStringObject struct is used to represent an SVG object, which consists of an object name and a type.
  *
  * The SVGStringObject struct contains the following members:
  * - object: A string representing the SVG element e.g. <rect x=1 ... />
  * - type: The type of object, used for sorting which object is displayed first.
  *
  * @see SVGType
  */
 struct SVGStringObject {
  SVGStringObject(const std::string &object, int priority): object(object), priority(priority) {
  }

  std::string object;
  int priority;
 };

 /**
  * @brief Represents a node in the SNFG tree structure.
  *
  * The SNFGNode class represents a node in the SNFG (Standard Nomenclature For Glycans) tree structure.
  * Each SNFGNode object can have a parent node, a sugar object, children nodes,
  * and various properties that define its position in the tree.
  * The class provides methods to navigate the tree structure and access the node's properties.
  */
 struct SNFGNode {
  SNFGNode() = default;

  explicit SNFGNode(Sugar *sugar) : sugar(sugar) {
  }

  /**
   @brief
   The color of an SNFG element
   */
  std::string snfg_colour;

  /**
  @brief
  The shape of an SNFG element
  */
  std::string snfg_shape;


  /**
  @brief Pointer to the linkage which describes the bond between the previous node and this node. Used to draw text
  on certain nodes
  */
  Linkage* linkage;


  /**
  @brief Pointer to the sugar which corresponds to this SNFGNode
  */
  Sugar *sugar;

  /**
  @brief Pointer to the residue which corresponds to this SNFGNode
  */
  gemmi::Residue *residue;

  /**
   @brief Pointer to the chain which corresponds to this SNFGNode
   */
  gemmi::Chain *chain;

  /**
   * @brief Pointer to the parent SNFGNode of the current node.
   *
   * This pointer represents the parent node of the current node and
   * is used to navigate the SNFG tree structure. It points to the parent
   * node object in the SNFG tree.
   *
   * @note The pointer may be nullptr if the current node is the root node
   *       or if it hasn't been assigned a parent yet.
   */
  SNFGNode *parent;
  /**
   * @brief A vector of unique pointers to SNFGNode objects representing the children nodes of the current SNFGNode.
   *
   * The children vector holds the child nodes of the current SNFGNode. Each child node is represented by a unique
   * pointer to an SNFGNode object.
   *
   * @note It is important to note that ownership of the child nodes is transferred to the parent node through
   * the use of unique pointers. This means that when a child node is added to the children vector, the parent
   * node becomes the exclusive owner of that child node. It is responsible for managing its lifetime and
   * eventually releasing the memory when it is no longer needed.
   *
   * @note By default, the children vector is empty, indicating that the current SNFGNode has no child nodes.
   * Child nodes can be added to the vector using the std::vector::push_back() function.
   */
  std::vector<std::unique_ptr<SNFGNode> > children;
  /**
   * @brief The position of the current node in its parent's children vector.
   *
   * The parent_position variable represents the position of the current node in its parent's children
   * vector. It is used to determine the node's relationship with its siblings (left, right, etc.).
   * The initial value of parent_position is 0, which indicates that the node is the leftmost child of its parent.
   *
   * @remark This variable is only valid if the parent of the node is correctly set and it has a valid parent position.
   */
  int parent_position = 0;
  /**
   * @brief y position of node
   */
  int y = 0;
  /**
   * @brief x position of node
   */
  int x = 0;
  /**
   * @brief y position modifier used to move entire subtreesz
   */
  int mod = 0;

  /**
   * @brief Checks if the current SNFGNode is the leftmost node among its siblings.
   *
   * This method checks if the current node is the leftmost node among its siblings
   * by comparing its position (parent_position) with 0.
   *
   * @return Returns true if the current SNFGNode is the leftmost node among its siblings, false otherwise.
   */
  [[nodiscard]] bool is_leftmost() const {
   return parent_position == 0;
  }

  /**
   * @brief Checks if the current SNFGNode is a leaf node (i.e., has no children).
   *
   * @return Returns true if the current SNFGNode is a leaf node, false otherwise.
   */
  [[nodiscard]] bool is_leaf() const {
   return children.empty();
  }

  /**
   * @brief Get the left sibling of the current node.
   *
   * If the node has no left sibling, nullptr is returned.
   *
   * @return The left sibling of the current node, or nullptr if there is no left sibling.
   */
  [[nodiscard]] SNFGNode *get_left_sibling() const {
   if (parent_position == 0) return nullptr;
   return parent->children[parent_position - 1].get();
  }

  /**
   * @brief Returns the right sibling of the current node.
   *
   * If the node has no left sibling, nullptr is returned.
   *
   * @return A pointer to the right sibling of the current node, or nullptr if there is no right sibling.
   */
  [[nodiscard]] SNFGNode *get_right_sibling() const {
   if (parent_position + 1 >= parent->children.size()) return nullptr;
   return parent->children[parent_position + 1].get();
  }

  /**
   * @brief Retrieve the leftmost sibling of the current node.
   *
   * This method returns a pointer to the leftmost sibling of the current node.
   * A sibling is a node that shares the same parent. If the current node does not
   * have any siblings or if the current node is the leftmost sibling, nullptr is returned.
   *
   * @return Pointer to the leftmost sibling of the current node, or nullptr if the current node has no leftmost sibling.
   *
   * @remark The returned pointer is owned by the parent node and should not be deleted by the caller.
   */
  [[nodiscard]] SNFGNode *get_leftmost_sibling() const {
   return parent->children[0].get();
  }

  /**
   * @brief Get the leftmost child of the current node. If there are no children, nullptr is returned.
   *
   * @return A pointer to the leftmost child SNFGNode object. If the current node has no children, nullptr is returned.
   * @remark The returned pointer is owned by the parent node and should not be deleted by the caller.
   */
  [[nodiscard]] SNFGNode *get_leftmost_child() const {
   if (this->children.empty()) return nullptr;
   return this->children[0].get();
  }

  /**
   * @brief Get the rightmost sibling node of the current node.
   *
   * @return The rightmost sibling node.
   * @remark The returned pointer is owned by the parent node and should not be deleted by the caller.
   */
  [[nodiscard]] SNFGNode *get_rightmost_sibling() const {
   return parent->children[parent->children.size() - 1].get();
  }

  /**
   * @brief Get the rightmost child of the current SNFGNode. If there are no children, nullptr is returned.
   *
   * @return A pointer to the rightmost child of the current SNFGNode.
   * @remark The returned pointer is owned by the parent node and should not be deleted by the caller.
   */
  [[nodiscard]] SNFGNode *get_rightmost_child() const {
   if (this->children.empty()) return nullptr;
   return this->children[this->children.size() - 1].get();
  }

  /**
   * @brief Checks if the current node has a right sibling.
   *
   * This method checks if the current node has a right sibling by comparing its position
   * in the parent's children vector with the size of the children vector minus one.
   *
   * @return True if the current node has a right sibling, false otherwise.
   */
  [[nodiscard]] bool has_right_sibling() const {
   return parent_position < parent->children.size() - 1;
  }

  /**
   * @brief Checks if the current node has a left sibling.
   *
   * @return True if the current node has a left sibling, False otherwise.
   *
   * @note This method assumes that the node's parent is correctly set and
   *       it has a valid parent position. If the parent position is 0, it
   *       indicates that the node does not have a left sibling.
   */
  [[nodiscard]] bool has_left_sibling() const {
   return parent_position != 0;
  }
 };

 // https://reingold.co/tidier-drawings.pdf
 /**
  * @class SNFG
  * @brief Class for creating SNFG (Symbol Nomenclature for Glycans) diagrams.
  *
  * The SNFG class is responsible for creating SNFG diagrams by taking a glycan structure,
  * a base residue, and using the Reingold-Tilford Positioning Algorithm to calculate the
  * positions of the nodes. Reference: https://reingold.co/tidier-drawings.pdf
  *
  * @note The SNFG class depends on the gemmi::Structure and ResidueDatabase classes.
  */
 class SNFG {
  /**
   * @brief Constructs a SNFG (Symbol Nomenclature for Glycans) object.
   *
   * The SNFG class is responsible for representing the symbolic nomenclature for glycans. It contains a constructor
   * that takes a gemmi::Structure pointer and a ResidueDatabase pointer as parameters.
   *
   * @param structure A pointer to a gemmi::Structure object that represents the structure of the glycans.
   * @param database A pointer to a ResidueDatabase object that contains the residue information for the glycans.
   */
 public:
  explicit SNFG(gemmi::Structure *structure, ResidueDatabase *database): m_structure(structure),
                                                                         m_database(database) {
  }

  /**
   * @brief Creates a Symbol Notation for Glycans (SNFG) representation of a glycan structure.
   *
   * The create_snfg method takes a Glycan object and a Glycosite object as parameters and generates a Symbol Notation for
   * Glycans (SNFG) representation of the glycan structure. It uses the provided glycan and base residue to create a hierarchical
   * tree structure of SNFGNodes representing the glycan structure. It then calculates the positions of the nodes in the tree,
   * prints the tree, creates SVG objects, orders the SVG objects, and finally generates an SVG file.
   *
   * @param glycan The Glycan object representing the glycan structure.
   * @param base_residue The Glycosite object representing the base residue of the glycan structure.
   *
   * @return The formatted SVG in string form.
   */
  std::string create_snfg(Glycan &glycan, Glycosite &base_residue);

 private:
  /**
   * @brief Calculates the positions of nodes within a SNFG (Symbol Nomeclature for Glycans) hierarchy.
   *
   * The calculate_node_positions method calculates the positions of nodes within a SNFG hierarchy. It takes in a pointer to
   * a SNFGNode object, representing the root node of the hierarchy.
   *
   * This method performs the following steps:
   * 1. Calls the initialise_nodes method to initialize the positions of all nodes in the hierarchy.
   * 2. Calls the calculate_initial_positions method to calculate the initial positions of all nodes.
   * 3. Calls the check_all_children_on_screen method to ensure all child nodes are within the screen bounds.
   * 4. Calls the calculate_final_positions method to calculate the final positions of all nodes.
   *
   * @param node A pointer to a SNFGNode object representing the root node of the SNFG hierarchy.
   *
   * @see initialise_nodes
   * @see calculate_initial_positions
   * @see check_all_children_on_screen
   * @see calculate_final_positions
   */
  void calculate_node_positions(SNFGNode *node);

  /**
   * @brief Initialises the nodes in the SNFG (Symbol Nomeclature for Glycans).
   *
   * The initialise_nodes method is used to initialise the nodes in a SNFG by setting their properties
   * such as the y-coordinate, x-coordinate, and mod value.
   *
   * This method takes two parameters:
   * - node: A pointer to the SNFGNode object that will be initialised.
   * - depth: An integer representing the depth level of the node in the SNFG.
   *
   * The method sets the y-coordinate of the node to -1, x-coordinate to the depth value, and mod value to 0.
   * It then recursively calls the initialise_nodes method for each child of the node.
   */
  void initialise_nodes(Sails::SNFGNode *node, int depth);

  /**
   * @brief Calculates the initial positions of the SNFG (Symbol Nomeclature for Glycans) nodes.
   *
   * The calculate_initial_positions method is used to calculate the initial positions of the SNFG nodes. It recursively
   * traverses the SNFG tree starting from the given node and assigns a y-coordinate to each node based on its position in
   * the tree.
   *
   * The y-coordinate of a node is determined by the following rules:
   * - If the node is a leaf node:
   *   - If it is not the leftmost sibling, its y-coordinate is set to the y-coordinate of its left sibling plus the size
   *     of the node and the distance between the siblings.
   *   - If it is the leftmost sibling, its y-coordinate is set to 0.
   * - If the node has only one child:
   *   - If it is the leftmost sibling, its y-coordinate is set to the y-coordinate of its leftmost child.
   *   - If it is not the leftmost sibling, its y-coordinate is set to the y-coordinate of its left sibling plus the size
   *     of the node and the distance between the siblings, and its mod (a variable used for optimizing node positioning)
   *     is set to the difference between its y-coordinate and the y-coordinate of its leftmost child.
   * - If the node has more than one child:
   *   - The y-coordinate of the node is set to the average of the y-coordinates of its leftmost and rightmost children.
   *   - If it is not the leftmost sibling, its y-coordinate is set to the y-coordinate of its left sibling plus the size
   *     of the node and the distance between the siblings, and its mod is set to the difference between its y-coordinate
   *     and the average of the y-coordinates of its leftmost and rightmost children.
   *
   * The method also checks for conflicts between sibling nodes and resolves them by adjusting their positions if necessary.
   * This step is performed only if the node has children and it is not the leftmost sibling.
   *
   * @param node The node for which to calculate the initial position.
   */
  void calculate_initial_positions(Sails::SNFGNode *node);

  /**
   * @brief Calculates the final positions of SNFG nodes.
   *
   * The calculate_final_positions method is used to calculate the final positions of SNFG (Symbol Nomeclature for Glycans)
   * nodes. It recursively updates the y-coordinate of each node based on the modification sum and the mod value
   * of the node. It also traverses the child nodes and calls the method recursively for each child.
   *
   * @param node A pointer to the SNFGNode whose final position is to be calculated.
   * @param mod_sum The cumulative modification sum. It is passed recursively to update the y-coordinate of each node.
   *
   * @note The method assumes that the node and mod_sum are valid inputs.
   */
  void calculate_final_positions(SNFGNode *node, float mod_sum);


  /**
   * @brief Checks for conflicts between a given node and its siblings in a tree.
   *
   * The check_for_conflicts method checks for conflicts between a given node and its siblings in an SNFG (Symbol Nomeclature for Glycans)
   * tree. It adjusts the position of the given node and its siblings if there is a conflict, ensuring that
   * nodes do not overlap.
   *
   * @param node A pointer to the SNFGNode object representing the node to check for conflicts.
   */
  void check_for_conflicts(SNFGNode *node);

  /**
   * @brief Centers the nodes between two given SNFGNodes.
   *
   * The center_nodes_between function is used to calculate the appropriate positions for the nodes between two given
   * SNFGNodes in an SNFG graph. It ensures that the nodes are evenly spaced between the two nodes and adjusts their
   * y-coordinates accordingly.
   *
   * @param lnode The left SNFGNode that defines the starting position of the nodes.
   * @param rnode The right SNFGNode that defines the ending position of the nodes.
   */
  void center_nodes_between(SNFGNode *lnode, SNFGNode *rnode);

  /**
   * @brief Checks and adjusts the position of all children nodes to ensure they are within the bounds of the screen.
   *
   * The check_all_children_on_screen method checks the position of all children nodes of the given SNFGNode to ensure
   * they are within the bounds of the screen. It adjusts the position of the nodes if necessary.
   *
   * The method performs the following steps:
   * 1. Calculates the left contour for all children nodes of the given node.
   * 2. Determines the amount of shift required for the nodes to be within the screen bounds by finding the minimum offset
   *    value among all nodes' right contour positions.
   * 3. Shifts the nodes vertically by the shift amount.
   *
   * @param node The SNFGNode for which to check and adjust the position of its children nodes.
   */
  void check_all_children_on_screen(SNFGNode *node);

  /**
   * @brief Calculates the left contour values for a given SNFGNode.
   *
   * The get_lcontour method calculates the left contour values for a given SNFGNode and updates the values parameter. The lcontour
   * values represent the minimum y-coordinate of each x-coordinate of the nodes in the SNFGNode tree. The mod_sum parameter
   * is used to calculate the modified sum of the y-coordinate values with the SNFGNode's mod value.
   *
   * @param node The SNFGNode for which to calculate the lcontour values.
   * @param mod_sum The modified sum of the y-coordinate values with the SNFGNode's mod value.
   * @param values A reference to a std::map<int, int> object to store the calculated lcontour values. The keys of the map
   *               represent the x-coordinates and the values represent the corresponding lcontour values.
   * @return None.
   *
   * @note The lcontour values are updated in the 'values' map parameter.
   * @note The 'mod_sum' parameter is modified within the function during recursion.
   */
  void get_lcontour(SNFGNode *node, int mod_sum, std::map<int, int> &values);

  /**
   * @brief Calculates the right contour values for each node in the SNFG.
   *
   * The get_rcontour method calculates and populates the right contour values for each node in the SNFG (Symbol Nomeclature for Glycans).
   * The rcontour value is the maximum y-coordinate value at each unique x-coordinate value in the SNFG.
   *
   * @param node The SNFGNode pointer representing the current node.
   * @param mod_sum The sum of the modifier values for the current node and its ancestors.
   *
   * @param values The reference to a std::map<int, int> that will store the calculated R contour values. The map is modified by
   * this method to add or update the contour values for each unique x-coordinate.
   * The key represents the unique x-coordinate value, and the value represents the maximum y-coordinate value for the given x-coordinate.
   *
   * @note The method modifies the values map by adding or updating the contour values for each unique x-coordinate.
   *
   * @note The method uses recursion to traverse the SNFG tree structure and calculate the R contour values for each node.
   */
  void get_rcontour(SNFGNode *node, int mod_sum, std::map<int, int> &values);


  /**
   * @brief Creates SVG (Scalable Vector Graphics) elements for the SNFG (Symbol Nomenclature for Glycans) nodes and their parent-child relationships.
   *
   * The create_svg method creates SVG elements based on the provided SNFG nodes and their parent-child relationships.
   * The SVG elements are stored in the snfg_objects vector.
   *
   * @param snfg_objects A reference to a vector of SVG objects. The SVG elements created by this method will be stored in this vector.
   * @param parent A pointer to the parent SNFGNode.
   * @param node A pointer to the current SNFGNode.
   *
   * @details The method calculates the coordinates of the current node and its parent, and creates SVG objects based on the node's SNFGShape and SNFGColour properties.
   * A line SVG object is also created using the create_svg_line method to represent the parent-child relationship between the current node and its parent.
   *
   * The method recursively calls itself for each child of the current node to create SVG elements for the entire SNFG structure.
   */
  void create_svg(std::vector<SVGStringObject> &snfg_objects, SNFGNode *parent, SNFGNode *node);

  /**
   * @brief Orders a vector of SVG objects by their type.
   *
   * The order_svg method takes a vector of SVG objects as input and orders them in descending order based on their type.
   * The type of an SVG object is determined by its "type" property.
   *
   * @param objects The vector of SVG objects to be ordered.
   */
  void order_svg(std::vector<Sails::SVGStringObject> &objects);

  /**
   * @brief Prints the tree structure starting from the given root node.
   *
   * The printTree method is used to print the tree structure starting from the specified root node. The tree is traversed recursively in a depth-first manner.
   *
   * @param root The root node of the tree to be printed.
   * @param node The current node being processed.
   * @param level The level of the current node in the tree.
   *
   * @note The tree printing is done in a depth-first manner, which means that the nodes are printed in the order they are visited.
   */
  void printTree(SNFGNode *root, SNFGNode *node, int level);

  void assign_snfg_types(Sugar *child,
                         SNFGNode *new_child);

  /**
   * @brief Forms the SNFG (Symbol Nomenclature for Glycans) node system.
   *
   * The form_snfg_node_system method forms the SNFG node system for a given root node, sugar, and glycans.
   *
   * @param root The root node of the SNFG node system.
   * @param sugar The sugar element used to form the node system.
   * @param glycan The glycans used to form the node system.
   */
  void form_snfg_node_system(SNFGNode *root, Sugar *sugar, Glycan &glycan);

  /**
   * @brief Creates the SVG header for an SVG document.
   *
   * This method generates the SVG header for an SVG (Scalable Vector Graphics) document by combining the SVG width and
   * height specified by the class constant values with the XML namespace declaration.
   *
   * The generated SVG header will be in the following format:
   *
   *      <svg width="[SVG_WIDTH]" height="[SVG_HEIGHT]" xmlns="http://www.w3.org/2000/svg">
   *
   * @note The SVG_WIDTH and SVG_HEIGHT are class constant values defined in the Sails::SNFG class.
   *
   * @return The SVG header as a std::string.
   */
  [[nodiscard]] std::string create_svg_header() const;

  /**
   * @brief Creates the footer of an SVG document.
   *
   * The create_svg_footer function generates the closing tag for the SVG element.
   *
   * The function returns a string representing the closing tag: "</svg>\n".
   *
   * @return String representing the closing tag for the SVG element.
   */
  [[nodiscard]] static std::string create_svg_footer();

  /**
   * @brief Creates an SVG line element.
   *
   * The create_svg_line method creates an SVG line element with the specified coordinates and returns it as an SVGObject.
   *
   * @param x1 - The x-coordinate of the starting point of the line.
   * @param y1 - The y-coordinate of the starting point of the line.
   * @param x2 - The x-coordinate of the ending point of the line.
   * @param y2 - The y-coordinate of the ending point of the line.
   *
   * @return An SVGObject representing the created SVG line element.
   *
   * Example usage:
   *
   * ```cpp
   * Sails::SVGObject line = Sails::SNFG::create_svg_line(10, 20, 30, 40);
   * ```
   */
  [[nodiscard]] static Sails::SVGStringObject create_svg_line(int x1, int y1, int x2, int y2);

  /**
   * @brief Creates an SVG text element with the given position and text.
   *
   * This function creates and returns an SVG text element with the specified position (x, y) and the provided text.
   * The created text element is defined by its position, font family, font size, fill color, and alignment.
   *
   * @param x The x-coordinate of the text element.
   * @param y The y-coordinate of the text element.
   * @param text The text content of the text element.
   * @return An SVGObject representing the created SVG text element.
   *
   * @note The font family, font size, fill color, and alignment are hardcoded in this implementation as "Verdana",
   *       20, black, and middle, respectively. Modify the implementation if different values are required.
   *
   * @see SVGObject
   * @see SVGType::text
   */
  [[nodiscard]] static Sails::SVGStringObject create_svg_text(int x, int y, const std::string &text);


  /**
   * @brief Create donor labels for a given SNFG node.
   *
   * The create_donor_labels method is used to create donor labels for a given SNFG (Symbol Nomenclature For Glycans) node.
   * Donor labels provide information about the linkage between glycan residues.
   *
   * This method takes three parameters: parent, node, and linkage.
   *
   * @param parent A pointer to the parent SNFG node.
   * @param node A pointer to the SNFG node for which donor labels need to be created.
   * @param linkage A pointer to the Linkage object that represents the linkage information.
   *
   * @note The parent parameter should not be null, and the node and linkage parameters can be null if there is no corresponding
   * linkage information available.
   *
   * @see SNFGNode, Linkage
   */
  [[nodiscard]] SVGStringObject create_donor_labels(SNFGNode * parent, SNFGNode * node, Linkage * linkage);


 private:
  /**
   * @brief Pointer to a gemmi::Structure object
   */
  gemmi::Structure *m_structure;

  /**
   * @brief Pointer to an instance of ResidueDatabase.
   */
  ResidueDatabase *m_database;

  /**
   * @brief The width of the SVG (Scalable Vector Graphics) document.
   */
  const int SVG_WIDTH = 1200;
  /**
   * @brief The height of the SVG (Scalable Vector Graphics) document.
   */
  const int SVG_HEIGHT = 800;

  /**
   * @brief Represents the size of a node.
   */
  const int node_size = 100;

  /**
   * @brief Represents the distance between sibling elements.
   */
  const int sibling_distance = 1;

  /**
   * @brief Represents the distance between two tree nodes.
   */
  const int tree_distance = 2;
 };
}
#endif //SAILS_SNFG_H
