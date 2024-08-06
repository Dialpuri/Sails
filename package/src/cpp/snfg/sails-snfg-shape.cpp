//
// Created by Jordan Dialpuri on 06/08/2024.
//

#include "../../include/snfg/sails-snfg-shape.h"

Sails::SVGStringObject Sails::SNFGShapeBase::draw() const {
    std::stringstream object;
    auto types = get_type();
    auto kwargs = get_kwargs();
    for (int i = 0; i < types.size(); i++) {
        object << "<" + types[i] + " ";
        for (auto& [k, v]: kwargs[i]) {
            object << k + "=\"" << v << "\" ";
        }
        object << "/>\n";
    }
    return {object.str(), get_priority()};
}

Sails::Circle::Circle(const SNFGNode *node) {
    kwargs["cx"] = std::to_string(node->x);
    kwargs["cy"] =  std::to_string(node->y) ;
    kwargs["r"]  = "30";
    kwargs["stroke"] = "black";
    kwargs["stroke-width"] = "4" ;
    kwargs["fill"] = node->snfg_colour;

    type="circle";

    priority = 1;
}


Sails::Square::Square(const SNFGNode *node) {
    kwargs["x"] = std::to_string(node->x-25);
    kwargs["y"] =  std::to_string(node->y-25) ;
    kwargs["width"]  = "50";
    kwargs["height"]  = "50";
    kwargs["stroke"] = "black";
    kwargs["stroke-width"] = "4" ;
    kwargs["fill"] = node->snfg_colour;

    type="rect";

    priority = 1;
}


Sails::Rectangle::Rectangle(const SNFGNode *node) {
    // border
    kwargs1["x"] = std::to_string(node->x-25);
    kwargs1["y"] =  std::to_string(node->y-25) ;
    kwargs1["width"]  = "160";
    kwargs1["height"]  = "50";
    kwargs1["stroke"] = "black";
    kwargs1["stroke-width"] = "4" ;
    kwargs1["fill"] = "white";
    kwargs1["rx"] = "10";
    kwargs1["ry"] = "10";
    type1="rect";

    // line
    kwargs2["x1"] = std::to_string(node->x+10);
    kwargs2["y1"] = std::to_string(node->y-25);
    kwargs2["x2"] = std::to_string(node->x+10);
    kwargs2["y2"] = std::to_string(node->y+25);
    kwargs2["stroke"] = "black";
    kwargs2["stroke-width"] = "4" ;
    type2="line";

    priority = 1;
}

std::unique_ptr<Sails::SNFGShapeBase> Sails::get_svg_shape(SNFGNode *node) {

    if (node->snfg_shape == "circle")
        return std::make_unique<Circle>(node);
    if (node->snfg_shape == "square")
        return std::make_unique<Square>(node);
    if (node->snfg_shape == "rectangle")
        return std::make_unique<Rectangle>(node);

}
