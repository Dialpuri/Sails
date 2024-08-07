//
// Created by Jordan Dialpuri on 07/08/2024.
//

#include "../../../include/snfg/shapes/square.h"

Sails::Square::Square(const SNFGNode *node) {
    kwargs["x"] = std::to_string(node->x - 25);
    kwargs["y"] = std::to_string(node->y - 25);
    kwargs["width"] = "50";
    kwargs["height"] = "50";
    kwargs["stroke"] = "black";
    kwargs["stroke-width"] = "4";
    kwargs["fill"] = node->snfg_colour;

    type = "rect";

    tooltip = format_tooltip(node);

    priority = 1;
}
