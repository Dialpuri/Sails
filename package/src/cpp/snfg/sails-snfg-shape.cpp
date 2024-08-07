//
// Created by Jordan Dialpuri on 06/08/2024.
//

#include "../../include/snfg/sails-snfg-shape.h"
#include "../../include/snfg/shapes/circle.h"
#include "../../include/snfg/shapes/rectangle.h"
#include "../../include/snfg/shapes/square.h"

std::string Sails::SNFGShapeBase::kwargs_to_string(std::map<std::string, std::string> kwargs) {
    std::stringstream object;
    for (auto &[k, v]: kwargs) {
        object << k + "=\"" << v << "\" ";
    }
    return object.str();
}

Sails::SVGStringObject Sails::SNFGShapeBase::draw() const {
    std::stringstream object;
    auto types = get_type();
    auto kwargs = get_kwargs();
    for (int i = 0; i < types.size(); i++) {
        object << "<" + types[i] + " ";
        object << kwargs_to_string(kwargs[i]);
        object << "/>\n";
    }
    object << get_special_tags();

    return {object.str(), get_priority()};
}


std::unique_ptr<Sails::SNFGShapeBase> Sails::get_svg_shape(SNFGNode *node) {
    if (node->snfg_shape == "circle")
        return std::make_unique<Sails::Circle>(node);
    if (node->snfg_shape == "square")
        return std::make_unique<Square>(node);
    if (node->snfg_shape == "rectangle")
        return std::make_unique<Rectangle>(node);
}
