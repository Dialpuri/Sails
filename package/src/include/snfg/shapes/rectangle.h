//
// Created by Jordan Dialpuri on 07/08/2024.
//

#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "../sails-snfg-shape.h"

namespace Sails {
    struct Rectangle final : SNFGShapeBase {
        static std::string get_donor_tag(const Sails::SNFGNode *node);

        static std::string get_donor_atom_tag(const SNFGNode *node);

        explicit Rectangle(const SNFGNode *node);

    protected:
        [[nodiscard]] KwargList get_kwargs() const override { return {rectangle_kwargs, line_kwargs}; }

        [[nodiscard]] std::vector<std::string> get_type() const override { return {type1, type2}; }

        [[nodiscard]] int get_priority() const override { return priority; }

        [[nodiscard]] std::string get_special_tags() const override {return special_obj;}

        [[nodiscard]] std::string get_tooltips() const override { return tooltip;}

    private:
        static std::string get_donor_string(const SNFGNode *node);


        std::map<std::string, std::string> rectangle_kwargs = {};
        std::map<std::string, std::string> line_kwargs = {};
        std::string special_obj;

        std::string type1;
        std::string type2;

        std::string tooltip;

        int priority;
    };
}

#endif //RECTANGLE_H
