//
// Created by Jordan Dialpuri on 04/08/2024.
//

#include "../include/sails-snfg.h"

std::string Sails::SNFG::create_svg_header() const {
    return "<svg width=\"" + std::to_string(SVG_WIDTH) + "\" height=\"" + std::to_string(SVG_HEIGHT) +
           "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
}

std::string Sails::SNFG::create_svg_footer() {
    return "</svg>\n";
}

std::string Sails::SNFG::create_svg_circle(int cx, int cy, int r, const std::string &color) {
    return "<circle cx=\"" + std::to_string(cx) + "\" cy=\"" + std::to_string(cy) + "\" r=\"" + std::to_string(r) +
           "\" stroke=\"black\" stroke-width=\"2\" fill=\"" + color + "\" />\n";
}

std::string Sails::SNFG::create_svg_line(int x1, int y1, int x2, int y2) {
    return "<line x1=\"" + std::to_string(x1) + "\" y1=\"" + std::to_string(y1) + "\" x2=\"" + std::to_string(x2) +
           "\" y2=\"" + std::to_string(y2) + "\" stroke=\"black\" stroke-width=\"2\" />\n";
}

std::string Sails::SNFG::create_svg_text(int x, int y, const std::string &text) {
    return "<text x=\"" + std::to_string(x) + "\" y=\"" + std::to_string(y) +
           "\" font-family=\"Verdana\" font-size=\"20\" fill=\"black\" text-anchor=\"middle\">" + text + "</text>\n";
}

void Sails::SNFG::printTree(SNFGNode* root, Sails::SNFGNode* node, int level) {
    if (node == nullptr) return;

    std::cout << "Level: " << level << " -> node: " << Sails::Utils::format_residue_from_site(node->sugar->site, m_structure) << ", x: " << node->x << ", y: " << node->y << std::endl;

    for (auto& child : node->children) {
        printTree(root, child.get(), level+1);
    }
}


void Sails::SNFG::form_snfg_node_system(SNFGNode *root, Sugar *sugar, Glycan &glycan) {
    int pos = 0;
    for (auto &child: glycan.adjacency_list[sugar]) {
        auto new_child = std::make_unique<SNFGNode>(child);
        new_child->parent_position = pos;
        new_child->parent = root;

        root->children.push_back(std::move(new_child));

        form_snfg_node_system(root->children.back().get(), child, glycan);
        pos++;
    }
}
void Sails::SNFG::calculate_node_positions(Sails::SNFGNode *node) {
    init_nodes(node, 0);

    calculate_initial_x(node);

    check_all_children_on_screen(node);

    calculate_final_positions(node, 0);
}

void Sails::SNFG::init_nodes(Sails::SNFGNode *node, int depth) {
    node->x = -1;
    node->y = depth;
    node->mod = 0;

    for (auto& child: node->children) {
        init_nodes(child.get(), depth+1);
    }
}

void Sails::SNFG::calculate_final_positions(Sails::SNFGNode *node, float mod_sum) {
    node->x += mod_sum;
    mod_sum += node->mod;

    for (auto& child: node->children) {
        calculate_final_positions(child.get(), mod_sum);
    }

    if (node->children.empty()) {

    }
}

void Sails::SNFG::calculate_initial_x(Sails::SNFGNode *node) {
    for (auto& child: node->children) {
        calculate_initial_x(child.get());
    }

    if (node->is_leaf()) {
        if (!node->is_leftmost()) {
            node->x = node->get_left_sibling()->x + node_size + sibling_distance;
        } else {
            node->x = 0;
        }
    } else if (node->children.size() == 1) {
        if (node->is_leftmost()) {
            node->x = node->get_leftmost_child()->x;
        } else {
            node->x = node->get_left_sibling()->x + node_size + sibling_distance;
            node->mod = node->x - node->get_leftmost_child()->x;
        }
    } else {
        int lc = node->get_leftmost_child()->x;
        int rc = node->get_rightmost_child()->x;

        std::cout << Utils::format_residue_from_site(node->sugar->site, m_structure) << " LC RC: " << lc << " " << rc << std::endl;
        int mid = (lc + rc) / 2;
        if (node->is_leftmost()) {
            node->x = mid;
        } else {
            node->x = node->get_left_sibling()->x + node_size + sibling_distance;
            node->mod = node->x - mid;
        }
    }
    if (node->children.size() > 0 && !node->is_leftmost()) {
        check_for_conflicts(node);
    }
}


void Sails::SNFG::check_for_conflicts(SNFGNode* node) {
    int min_distance = tree_distance + node_size;
    int shift = 0;

    std::map<int, int> node_contour;
    get_lcontour(node, 0, node_contour);

    auto sibling = node->get_leftmost_sibling();
    while(sibling != nullptr && sibling != node) {
        std::map<int, int> sibling_contour;
        get_rcontour(sibling, 0, sibling_contour);

        auto last_sibling_contour = std::prev(sibling_contour.end())->first;
        auto last_node_contour = std::prev(node_contour.end())->first;
        for (int i = node->y + 1; i <= std::min(last_sibling_contour, last_node_contour); i++) {
            int distance = node_contour[i] - sibling_contour[i];
            if (distance + shift < min_distance) {
                shift = min_distance - distance;
            }
        }

        if (shift > 0) {
            node->x += shift;
            node->mod += shift;

            center_nodes_between(node, sibling);
            shift = 0;
        }
        sibling = sibling->get_right_sibling();
    }
}

void Sails::SNFG::center_nodes_between(Sails::SNFGNode *lnode, Sails::SNFGNode *rnode) {
    int li = lnode->parent_position;
    int ri = rnode->parent_position;

    int node_between = (ri - li) - 1;
    if (node_between > 0) {
        int node_distance = (lnode->x - rnode->x) / (node_between+1);
        int count = 1;
        for (int i = li+1; i < ri; i++) {
            SNFGNode* mid_node = lnode->parent->children[i].get();
            int desired_x = rnode->x + (node_distance*count);
            int offset = desired_x - mid_node->x;
            mid_node->x += offset;
            mid_node->x += offset;
            count++;
        }
        check_for_conflicts(lnode);
    }
}

void Sails::SNFG::check_all_children_on_screen(Sails::SNFGNode *node) {
    std::map<int, int> node_contour;

    get_lcontour(node, 0, node_contour);

    int shift_amount = 0;
    for (auto& [k, v]: node_contour) {
        if (node_contour[k] + shift_amount < 0) {
            shift_amount = -node_contour[k];
        }
    }
    if (shift_amount > 0) {
        node->x += shift_amount;
        node->mod += shift_amount;
    }
}

void Sails::SNFG::get_lcontour(Sails::SNFGNode *node, int mod_sum, std::map<int, int> &values) {
    if (values.find(node->y) == values.end()) {
        values[node->y] = node->x + mod_sum;
    } else {
        values[node->y] = std::min(values[node->y], node->x + mod_sum);
    }

    mod_sum += node->mod;

    for (auto& child: node->children) {
        get_lcontour(child.get(), mod_sum, values);
    }
}

void Sails::SNFG::get_rcontour(Sails::SNFGNode *node, int mod_sum, std::map<int, int> &values) {
    if (values.find(node->y) == values.end()) {
        values[node->y] = node->x + mod_sum;
    } else {
        values[node->y] = std::max(values[node->y], node->x + mod_sum);
    }

    mod_sum += node->mod;

    for (auto& child: node->children) {
        get_rcontour(child.get(), mod_sum, values);
    }
}

void Sails::SNFG::create_svg(std::ofstream &f, SNFGNode *parent, SNFGNode *node) {
    int x = 400+node->x;
    int y = 200+100*node->y;
    int px = 400+parent->x;
    int py = 200+100*parent->y;
    f << create_svg_circle(x, y, 30, "lightblue");
    f << create_svg_text(x, y, Utils::format_residue_from_site(node->sugar->site, m_structure));
    f << create_svg_line(px, py, x, y);
    for (auto& child : node->children) {
        create_svg(f, node, child.get());
    }
}

std::string Sails::SNFG::create_snfg(Glycan &glycan, Glycosite &base_residue) {
    auto base_sugar = glycan.sugars[base_residue].get();

    auto root = std::make_unique<SNFGNode>(base_sugar);
    form_snfg_node_system(root.get(), base_sugar, glycan);

    calculate_node_positions(root.get());

    printTree(root.get(), root.get(), 0);


    std::ofstream f("tree.svg");
    f << create_svg_header();
    create_svg(f, root.get(), root.get());
    f << create_svg_footer();
    f.close();
    return "";
}

