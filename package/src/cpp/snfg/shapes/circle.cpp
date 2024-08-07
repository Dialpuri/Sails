//
// Created by Jordan Dialpuri on 07/08/2024.
//

#include "../../../include/snfg/shapes/circle.h"

Sails::Circle::Circle(const SNFGNode *node) {
    kwargs["cx"] = std::to_string(node->x);
    kwargs["cy"] = std::to_string(node->y);
    kwargs["r"] = "30";
    kwargs["stroke"] = "black";
    kwargs["stroke-width"] = "4";
    kwargs["fill"] = node->snfg_colour;

    type = "circle";

    tooltip = format_tooltip(node);

    priority = 1;
}
