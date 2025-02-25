//
// Created by Jordan Dialpuri on 07/08/2024.
//

#ifndef DIAMOND_H
#define DIAMOND_H

#include "../sails-snfg-shape.h"

namespace Sails {
    struct Diamond : SNFGShapeBase {
        explicit Diamond(const SNFGNode *node);

    protected:
        [[nodiscard]] KwargList get_kwargs() const override { return {kwargs}; };

        [[nodiscard]] std::vector<std::string> get_type() const override { return {type}; }

        [[nodiscard]] int get_priority() const override { return priority; }

        [[nodiscard]] std::string get_special_tags() const override { return ""; }; // no special tags for diamond

        [[nodiscard]] std::string get_tooltips() const override { return tooltip;}

    private:
        std::map<std::string, std::string> kwargs = {};
        std::string type;

        std::string tooltip;

        int priority;
    };
}
#endif //DIAMOND_H
