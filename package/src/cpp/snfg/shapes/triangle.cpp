//
// Created by Jordan Dialpuri on 07/08/2024.
//

#include "../../../include/snfg/shapes/triangle.h"

Sails::Triangle::Triangle(const SNFGNode *node) {
    kwargs["rx"] = std::to_string(node->x - 25);
    kwargs["ry"] = std::to_string(node->y - 25);

    float point1_x = node->x - 25;
    float point1_y = node->y + 25;

    float point2_x = node->x ;
    float point2_y = node->y - 25;

    float point3_x = node->x + 25;
    float point3_y = node->y + 25;

    std::string points;
    points += std::to_string(point1_x) + "," + std::to_string(point1_y) + " ";
    points += std::to_string(point2_x) + "," + std::to_string(point2_y) + " ";
    points += std::to_string(point3_x) + "," + std::to_string(point3_y);

    kwargs["points"] = points;
    kwargs["stroke"] = "black";
    kwargs["stroke-width"] = "4";
    kwargs["fill"] = node->snfg_colour;

    type = "polygon";

    tooltip = format_tooltip(node);

    priority = 1;
}
