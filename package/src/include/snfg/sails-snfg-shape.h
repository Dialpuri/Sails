//
// Created by Jordan Dialpuri on 06/08/2024.
//

#ifndef SAILS_SNFG_SHAPE_H
#define SAILS_SNFG_SHAPE_H
#include "sails-snfg.h"

namespace Sails {
    typedef std::vector<std::map<std::string, std::string>> KwargList;

    struct SNFGShapeBase {
        virtual ~SNFGShapeBase() = default;

        [[nodiscard]] SVGStringObject draw() const;

    protected:
        [[nodiscard]] virtual KwargList get_kwargs() const = 0;

        [[nodiscard]] virtual std::vector<std::string> get_type() const = 0;

        [[nodiscard]] virtual int get_priority() const = 0;
    };

    struct Circle final : SNFGShapeBase {
        explicit Circle(const SNFGNode *node);

    protected:
        [[nodiscard]] KwargList get_kwargs() const override { return {kwargs}; };

        [[nodiscard]] std::vector<std::string> get_type() const override { return {type}; }

        [[nodiscard]] int get_priority() const override { return priority; }


    private:
        std::map<std::string, std::string> kwargs = {};
        std::string type;
        int priority;
    };

    struct Square final : SNFGShapeBase {
        explicit Square(const SNFGNode *node);

    protected:
        [[nodiscard]] KwargList get_kwargs() const override { return {kwargs}; };

        [[nodiscard]] std::vector<std::string> get_type() const override { return {type}; }

        [[nodiscard]] int get_priority() const override { return priority; }

    private:
        std::map<std::string, std::string> kwargs = {};
        std::string type;
        int priority;
    };

    struct Rectangle final : SNFGShapeBase {
        explicit Rectangle(const SNFGNode *node);

    protected:
        [[nodiscard]] KwargList get_kwargs() const override { return {kwargs1, kwargs2}; };

        [[nodiscard]] std::vector<std::string> get_type() const override { return {type1, type2}; }

        [[nodiscard]] int get_priority() const override { return priority; }

    private:
        std::map<std::string, std::string> kwargs1 = {};
        std::map<std::string, std::string> kwargs2 = {};

        std::string type1;
        std::string type2;

        int priority;
    };


    std::unique_ptr<SNFGShapeBase> get_svg_shape(SNFGNode *node);
}


#endif //SAILS_SNFG_SHAPE_H
