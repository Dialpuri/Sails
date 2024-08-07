//
// Created by Jordan Dialpuri on 07/08/2024.
//

#include "../../../include/snfg/shapes/rectangle.h"

std::string Sails::Rectangle::get_donor_tag(const SNFGNode *node) {
    std::stringstream stream;
    std::map<std::string, std::string> donor_residue_name_kwargs;
    donor_residue_name_kwargs["x"] = std::to_string(node->x + 75);
    donor_residue_name_kwargs["y"] = std::to_string(node->y + 5);
    donor_residue_name_kwargs["fill"] = "black";
    donor_residue_name_kwargs["text-anchor"] = "middle";
    donor_residue_name_kwargs["font-weight"] = "bold";
    donor_residue_name_kwargs["font-family"] = "Helvetica";
    donor_residue_name_kwargs["font-size"] = "24";

    std::map<std::string, std::string> donor_seqid_kwargs;
    donor_seqid_kwargs["baseline-shift"] = "sub";
    donor_seqid_kwargs["font-weight"] = "normal";
    donor_seqid_kwargs["font-size"] = "20";

    stream << "<text ";
    stream << kwargs_to_string(donor_residue_name_kwargs);
    stream << ">";
    stream << node->residue->name;
    stream << "<tspan ";
    stream << kwargs_to_string(donor_seqid_kwargs);
    stream << ">";
    stream << node->chain->name << "/" << node->residue->seqid.num.str();
    stream << "</tspan></text>\n";

    return stream.str();
}

std::string Sails::Rectangle::get_donor_atom_tag(const SNFGNode *node) {
    std::map<std::string, std::string> donor_atom_kwargs;
    donor_atom_kwargs["x"] = std::to_string(node->x - 15);
    donor_atom_kwargs["y"] = std::to_string(node->y + 10);
    donor_atom_kwargs["font-weight"] = "bold";
    donor_atom_kwargs["font-family"] = "Helvetica";
    donor_atom_kwargs["font-size"] = "24";
    donor_atom_kwargs["stroke"] = "blue";

    std::stringstream stream;
    stream << "<text ";
    stream << kwargs_to_string(donor_atom_kwargs);
    stream << ">";
    stream << get_donor_string(node);
    stream << "</text>\n";
    return stream.str();
}

Sails::Rectangle::Rectangle(const SNFGNode *node) {
    // border
    rectangle_kwargs["x"] = std::to_string(node->x - 25);
    rectangle_kwargs["y"] = std::to_string(node->y - 25);
    rectangle_kwargs["width"] = "160";
    rectangle_kwargs["height"] = "50";
    rectangle_kwargs["stroke"] = "black";
    rectangle_kwargs["stroke-width"] = "4";
    rectangle_kwargs["fill"] = "white";
    rectangle_kwargs["rx"] = "10";
    rectangle_kwargs["ry"] = "10";
    type1 = "rect";

    // line
    line_kwargs["x1"] = std::to_string(node->x + 10);
    line_kwargs["y1"] = std::to_string(node->y - 25);
    line_kwargs["x2"] = std::to_string(node->x + 10);
    line_kwargs["y2"] = std::to_string(node->y + 25);
    line_kwargs["stroke"] = "black";
    line_kwargs["stroke-width"] = "4";
    type2 = "line";

    // N or C
    std::stringstream special_stream;
    special_stream << get_donor_atom_tag(node);
    special_stream << get_donor_tag(node);
    special_obj = special_stream.str();

    tooltip = format_tooltip(node);

    priority = 1;
}

std::string Sails::Rectangle::get_donor_string(const SNFGNode *node) {
    if (node->sugar == nullptr) { throw std::runtime_error("Attempting to draw rectangle with a null sugar"); }
    const std::string donor_atom = node->sugar->atom;
    return {donor_atom[0]};
}
