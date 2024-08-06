//
// Created by Jordan Dialpuri on 04/08/2024.
//

#ifndef SAILS_SNFG_H
#define SAILS_SNFG_H

#include <gemmi/model.hpp>

#include "sails-model.h"
#include "sails-telemetry.h"
#include "sails-topology.h"

namespace Sails {

    // struct SNFGNode {
    //     SNFGNode() = default;;
    //
    //     explicit SNFGNode(Sugar* sugar): sugar(sugar) {}
    //
    //     SNFGNode(Sugar *sugar, std::vector<std::unique_ptr<SNFGNode>> &children, int x, int y)
    //         : sugar(sugar),
    //           children(std::move(children)),
    //           x(x),
    //           y(y) {
    //     }
    //
    //     Sugar* sugar{};
    //     std::vector<std::unique_ptr<SNFGNode>> children;
    //     int x{};
    //     int y{};
    // };

    struct SNFGNode {
        SNFGNode() = default;

        explicit SNFGNode(Sugar* sugar) : sugar(sugar){}

        Sugar* sugar;
        SNFGNode* parent;
        std::vector<std::unique_ptr<SNFGNode>> children;
        int parent_position = 0;
        int x = 0;
        int y = 0;
        int mod = 0;
        int prelim_x = 0;
        int preLim_y = 0;

        [[nodiscard]] bool is_leftmost() const {
            return parent_position == 0;
        }

        [[nodiscard]] bool is_leaf() const {
            return children.empty();
        }

        [[nodiscard]] SNFGNode* get_left_sibling() const {
            if (parent_position == 0) return nullptr;
            return parent->children[parent_position-1].get();
        }

        [[nodiscard]] SNFGNode* get_right_sibling() const {
            if (parent_position+1 >= parent->children.size()) return nullptr;
            return parent->children[parent_position+1].get();
        }

        [[nodiscard]] SNFGNode* get_leftmost_sibling() const {
            return parent->children[0].get();
        }

        [[nodiscard]] SNFGNode* get_leftmost_child() const {
            return this->children[0].get();
        }

        [[nodiscard]] SNFGNode* get_rightmost_sibling() const {
            return parent->children[parent->children.size()-1].get();
        }

        [[nodiscard]] SNFGNode* get_rightmost_child() const {
            return this->children[parent->children.size()-1].get();
        }

        [[nodiscard]] bool has_right_sibling() const {
            return parent_position < parent->children.size()-1;
        }

        [[nodiscard]] bool has_left_sibling() const {
            return parent_position != 0;
        }

    };


    // Walker Algorithm - https://www.cs.unc.edu/techreports/89-034.pdf
    class SNFG {
    public:

        explicit SNFG(gemmi::Structure* structure): m_structure(structure) {}

        std::string create_snfg(Glycan& glycan, Glycosite& base_residue);

    private:

        // walker algorith,
        void first_walk(SNFGNode *root, SNFGNode* node, int depth);

        bool second_walk(SNFGNode *node, int level = 0, int modsum = 0);

        void create_svg(std::ofstream& f, SNFGNode* node);

        void printTree(SNFGNode *root, SNFGNode* node, int level);

        void form_snfg_node_system(SNFGNode* root, Sugar* sugar, Glycan& glycan);

        void check_for_conflicts(SNFGNode* node);

        void contour(SNFGNode* node, int mod_sum, std::map<int, int>& values);

        void center_nodes(SNFGNode* l, SNFGNode* r);

        void apportion(SNFGNode* root, SNFGNode* node, int depth);

        SNFGNode* get_leftmost_node(SNFGNode* node, int level, int depth);

        [[nodiscard]] SNFGNode* get_previous_node_at_level(SNFGNode *root, SNFGNode *node) const;

        [[nodiscard]] std::string create_svg_header() const ;

        [[nodiscard]] static std::string create_svg_footer() ;

        [[nodiscard]] static std::string create_svg_circle(int cx, int cy, int r, const std::string &color) ;

        [[nodiscard]] static std::string create_svg_line(int x1, int y1, int x2, int y2) ;

        [[nodiscard]] static std::string create_svg_text(int x, int y, const std::string &text) ;

    private:
        gemmi::Structure* m_structure;

        // Constants for SVG dimensions
        const int SVG_WIDTH = 1200;
        const int SVG_HEIGHT = 1200;
        const int NODE_RADIUS = 20;
        const int VERTICAL_SPACING = 80;
        const int HORIZONTAL_SPACING = 60;
    };


}
#endif //SAILS_SNFG_H
