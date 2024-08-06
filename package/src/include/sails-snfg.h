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


    enum SVGType {
        square, circle, line, text
    };

    struct SVGObject {
        SVGObject(const std::string& object, SVGType type): object(object), type(type) {}
        std::string object;
        SVGType type;
    };

    struct SNFGNode {
        SNFGNode() = default;

        explicit SNFGNode(Sugar* sugar) : sugar(sugar){}

        Sugar* sugar;
        SNFGNode* parent;
        std::vector<std::unique_ptr<SNFGNode>> children;
        int parent_position = 0;
        int y = 0;
        int x = 0;
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
            return this->children[this->children.size()-1].get();
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

        explicit SNFG(gemmi::Structure* structure, ResidueDatabase* database): m_structure(structure), m_database(database) {}

        std::string create_snfg(Glycan& glycan, Glycosite& base_residue);

    private:
        // reingold algorith,

        void calculate_node_positions(SNFGNode* node);

        void initialise_nodes(Sails::SNFGNode *node, int depth);

        void calculate_final_positions(SNFGNode* node, float mod_sum);

        void calculate_initial_positions(Sails::SNFGNode *node);

        void check_for_conflicts(SNFGNode* node);

        void center_nodes_between(SNFGNode* lnode, SNFGNode* rnode);

        void check_all_children_on_screen(SNFGNode* node);

        void get_lcontour(SNFGNode* node, int mod_sum, std::map<int, int>& values);

        void get_rcontour(SNFGNode* node, int mod_sum, std::map<int, int>& values);


        void create_svg(std::vector<Sails::SVGObject>& snfg_objects, SNFGNode *parent, SNFGNode *node);

        void order_svg(std::vector<Sails::SVGObject> &objects);


        void printTree(SNFGNode *root, SNFGNode* node, int level);

        void form_snfg_node_system(SNFGNode* root, Sugar* sugar, Glycan& glycan);


        [[nodiscard]] std::string create_svg_header() const ;

        [[nodiscard]] static std::string create_svg_footer() ;

        [[nodiscard]] static Sails::SVGObject create_svg_circle(int cx, int cy, int r, const std::string &color) ;

        [[nodiscard]] static Sails::SVGObject create_svg_square(int x, int y, int s, const std::string &color) ;

        [[nodiscard]] static Sails::SVGObject create_svg_line(int x1, int y1, int x2, int y2) ;

        [[nodiscard]] static Sails::SVGObject create_svg_text(int x, int y, const std::string &text) ;

    private:
        gemmi::Structure* m_structure;
        ResidueDatabase* m_database;

        // Constants for SVG dimensions
        const int SVG_WIDTH = 2400;
        const int SVG_HEIGHT = 2400;

        const int node_size = 100;
        const int sibling_distance = 1;
        const int tree_distance = 2;

    };


}
#endif //SAILS_SNFG_H
