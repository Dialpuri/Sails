//
// Created by Jordan Dialpuri on 06/08/2024.
//

#ifndef SAILS_SNFG_SHAPE_H
#define SAILS_SNFG_SHAPE_H

#include "sails-snfg.h"

namespace Sails {
 typedef std::vector<std::map<std::string, std::string> > KwargList;

 /**
  * @brief Base class for SNFG shape objects.
  *
  * This class provides a base for creating SNFG shape objects. It defines common functions and properties that derived classes must implement.
  * The draw() method creates the SVG object and produces elements based on the KwargList provided and then appends
  * any special tags (additional objects that are required to come along witth this SNFGShape).
  */
 struct SNFGShapeBase {
  /**
   * @brief Destructor for the SNFGShapeBase class.
   *
   * This destructor is declared as virtual and has a default implementation.
   * It is responsible for cleaning up any resources used by the SNFGShapeBase object.
   *
   * @note This destructor should be called when an instance of SNFGShapeBase is no longer needed.
   *
   * @see SNFGShapeBase
   */
  virtual ~SNFGShapeBase() = default;

  /**
   * @brief Converts key-value pairs from a map to a string representation.
   *
   * This method takes a map of key-value pairs and converts them to a string
   * representation. Each pair is formatted as "key=\"value\" " (note the space
   * at the end).
   *
   * @param kwargs The map containing the key-value pairs.
   * @return The string representation of the key-value pairs.
   */
  static std::string kwargs_to_string(std::map<std::string, std::string> kwargs);

  /**
   * @brief Retrieves the tooltips for a given SNFGNode.
   *
   * @return The tooltips for the node as a string.
   */
  static std::string format_tooltip(const SNFGNode *node);

  /**
   * @brief Draws the shape and returns it as an SVGStringObject.
   *
   * This method constructs an SVG string object representing the shape
   * according to the shape's type and keyword arguments. It then appends any special tags
   * specific to the shape. The resulting SVG string object is returned along with the
   * shape's priority as an SVGStringObject.
   *
   * @return The SVGStringObject representing the shape and its priority.
   */
  [[nodiscard]] SVGStringObject draw() const;

 protected:
  /**
   * @brief Get the keyword arguments.
   *
   * This method returns the keyword arguments as a KwargList object.
   *
   * @return The keyword arguments as a KwargList object.
   */
  [[nodiscard]] virtual KwargList get_kwargs() const = 0;

  /**
   * @brief Retrieves the type of the object.
   *
   * This method returns the type of the object as a vector of strings. The actual contents
   * of the vector may vary depending on the implementation.
   *
   * @return A vector of strings representing the type of the object.
   *
   * @note This method is virtual, and it must be implemented by all derived classes.
   *
   * @attention The returned vector may be empty if the object has no specific type.
   *
   * @see ObjectType
   */
  [[nodiscard]] virtual std::vector<std::string> get_type() const = 0;

  /**
   * @brief Get the priority of the object.
   *
   * This method is used to retrieve the priority of the object.
   *
   * @return The priority value of the object.
   */
  [[nodiscard]] virtual int get_priority() const = 0;

  /**
   * @brief Get the special tags.
   *
   * This method returns a string containing the special tags associated with the object.
   *
   * @return A string containing the special tags.
   */
  [[nodiscard]] virtual std::string get_special_tags() const = 0;


  /**
  * @brief Get the tool.
  *
  * This method returns a string containing the tooltips associated with the object.
  *
  * @return A string containing the tooltips.
  */
  [[nodiscard]] virtual std::string get_tooltips() const = 0;
 };


 /**
  * @brief Retrieves an SVG shape for the given SNFGNode.
  *
  * This method returns a unique_ptr to an SNFGShapeBase object, representing an SVG shape,
  * based on the shape specified in the SNFGNode parameter.
  *
  * @param node The SNFGNode for which to retrieve the SVG shape.
  *
  * @return A unique_ptr to an SNFGShapeBase object representing the SVG shape,
  *         or nullptr if the shape is not supported.
  */
 std::unique_ptr<SNFGShapeBase> get_svg_shape(SNFGNode *node);
}


#endif //SAILS_SNFG_SHAPE_H
