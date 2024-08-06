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

    auto s = get_previous_node_at_level(root, node);

    if (s != nullptr)
    std::cout << "leftmost node is " << Utils::format_residue_from_site(s->sugar->site, m_structure) << std::endl;

    for (auto& child : node->children) {
        printTree(root, child.get(), level+1);
    }
}

void Sails::SNFG::create_svg(std::ofstream& f, SNFGNode* node) {
    f << create_svg_circle(400+100*node->x, 200+100*node->y, 30, "lightblue");
    f << create_svg_text(400+100*node->x, 200+100*node->y, Utils::format_residue_from_site(node->sugar->site, m_structure));
    for (auto& child : node->children) {
        create_svg(f, child.get());
    }
}





void Sails::SNFG::check_for_conflicts(SNFGNode* node) {
    int min_distance = 0 ;
    int shift = 0;

    std::map<int, int> node_contour;
    contour(node, 0, node_contour);

    auto sibling = node->get_left_sibling();
    while(sibling != nullptr && sibling != node) {
        std::map<int, int> sibling_contour;
        contour(sibling, 0, sibling_contour);

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

            center_nodes(node, sibling);
            shift = 0;
        }
        sibling = sibling->get_right_sibling();
    }
}

void Sails::SNFG::contour(SNFGNode *node, int mod_sum, std::map<int, int> &values) {
    if (values.find(node->y) == values.end()) {
        values[node->y] = node->x+mod_sum;
    } else {
        values[node->y] = std::min(values[node->y], node->x+mod_sum);
    }
    mod_sum += node->mod;

    for (auto& child: node->children) {
        contour(child.get(), mod_sum, values);
    }
}

void Sails::SNFG::center_nodes(SNFGNode *l, SNFGNode *r) {

    int li = l->parent_position;
    int ri = r->parent_position;

    int node_between = (ri - li) - 1;
    if (node_between > 0) {
        int node_distance = (l->x - r->x) / (node_between+1);
        int count = 1;
        for (int i = li+1; i < ri; i++) {
            SNFGNode* mid_node = l->children[i].get();
            int desired_x = r->x + (node_distance*count);
            int offset = desired_x - mid_node->x;
            mid_node->x += offset;
            mid_node->x += offset;
            count++;
        }
        check_for_conflicts(l);
    }

}


std::string Sails::SNFG::create_snfg(Glycan &glycan, Glycosite &base_residue) {
    auto base_sugar = glycan.sugars[base_residue].get();

    auto root = std::make_unique<SNFGNode>(base_sugar);
    form_snfg_node_system(root.get(), base_sugar, glycan);

    first_walk(root.get(), root.get(), 0);

    // second_walk(root.get(), 0, 0);
    printTree(root.get(), root.get(), 0);

    std::ofstream f("tree.svg");
    f << create_svg_header();
    create_svg(f, root.get());
    f << create_svg_footer();
    f.close();
    return "";
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

Sails::SNFGNode * Sails::SNFG::get_previous_node_at_level(SNFGNode *root, SNFGNode *node) const {
    if (root == nullptr) throw std::runtime_error("NPR passed to root");

    std::queue<SNFGNode*> q;
    q.push(root);

    while (!q.empty()) {
        size_t level_width = q.size();
        SNFGNode* previous_node = nullptr;

        for (int i = 0; i < level_width; i++) {
            SNFGNode* current_node = q.front();
            q.pop();

            if (current_node == node) {
                return previous_node;
            }

            previous_node = current_node;

            for (auto& child: current_node->children) {
                q.push(child.get());
            }
        }
    }

    return nullptr;
}


void Sails::SNFG::first_walk(SNFGNode *root, SNFGNode *node, int depth) {

    auto left_neighbour = get_previous_node_at_level(root, node);
    node->mod = 0;

    if (node->is_leaf()) {
        if (node->get_left_sibling() == nullptr) {
            node->prelim_x = 0;
        } else {
            node->prelim_x = node->get_left_sibling()->prelim_x + 1;
        }
    } else {
        auto l = node->get_leftmost_child();
        auto r = l;

        while (r->has_right_sibling()) {
            r = r->get_right_sibling(); // check
            first_walk(root, r, depth + 1);
        }

        auto midpoint = l->prelim_x + r->prelim_x / 2;

        if (node->has_left_sibling()) {
            node->prelim_x = node->get_left_sibling()->prelim_x + 1;
            node->mod = node->prelim_x - midpoint;
            apportion(root, node, depth);
        }
        else {
            node->prelim_x = midpoint;
        }
    }

}


bool Sails::SNFG::second_walk(SNFGNode *node, int level, int modsum) {
    if (level < 100) {
        int xtmp = node->prelim_x + modsum;
        int ytmp = level;

        node->x = xtmp;
        node->y = ytmp;

        auto result = true;
        if (!node->children.empty()) {
            result = second_walk(node->get_leftmost_child(), level + 1, modsum + node->mod);
        }

        if (result && node->has_right_sibling()) {
            second_walk(node->get_right_sibling(), level + 1, modsum);
        }
        return false;
    }
    return true;
}


Sails::SNFGNode* Sails::SNFG::get_leftmost_node(SNFGNode* node, int level, int depth) {
    if (level >= depth) return node;
    if (node->is_leaf()) return nullptr;

    auto r = node->get_leftmost_child();
    auto l = get_leftmost_node(r, level + 1, depth);

    while (l == nullptr && r->has_right_sibling()) {
        r = r->get_right_sibling();
        l = get_leftmost_node(r, level+1, depth);
    }
    return l;
}


void Sails::SNFG::apportion(SNFGNode *root, SNFGNode *node, int depth) {
    auto l = node->get_leftmost_child();
    auto neighbour = get_previous_node_at_level(root, node);
    int compare_depth = 1;
    int stop_depth = 100-depth;

    while (l != nullptr && neighbour != nullptr && compare_depth <= stop_depth) {
        int l_modsum = 0;
        int r_modsum = 0;
        auto ancestor_leftmost = l;
        auto ancestor_neighbour = neighbour;

        for (int i = 0; i < compare_depth; i++) {
            ancestor_leftmost = ancestor_leftmost->parent;
            ancestor_neighbour = ancestor_neighbour->parent;
            r_modsum += ancestor_leftmost->mod;
            l_modsum += ancestor_neighbour->mod;
        }

        int move_distance = neighbour->prelim_x + l_modsum + 1 - l->prelim_x + r_modsum;
        if (move_distance > 0) {
            auto tmp = node;
            int left_siblings = 0;
            while (tmp != nullptr && tmp != ancestor_neighbour) {
                left_siblings++;
                tmp = node->get_left_sibling();
            }

            if (tmp != nullptr) {
                int portion = move_distance / left_siblings;
                tmp = node;
                while (tmp == ancestor_neighbour) {
                    tmp->prelim_x += move_distance;
                    tmp->mod += move_distance;
                    move_distance -= portion;
                    tmp = tmp->get_left_sibling();
                }
            }
        }

        compare_depth++;

        if (l->is_leaf()) {
            l = get_leftmost_node(node, 0, compare_depth);
        }
        else {
            l = l->get_leftmost_child();
        }
    }
}



// void Sails::SNFG::first_walk(SNFGNode *node, int depth) {
    // if (node == nullptr) return;
    //
    // for (auto& child: node->children)
    //     first_walk(child.get(), depth + 1);
    //
    // if (node->is_leaf()) {
    //     if (!node->is_leftmost()) {
    //         node->x = node->get_left_sibling()->x + 1;
    //     } else {
    //         node->x = 0;
    //     }
    // }
    // else if (node->children.size() == 1) {
    //     if (node->is_leftmost()) {
    //         node->x = node->children[0]->x;
    //     } else {
    //         node->x = node->get_left_sibling()->x + 1;
    //         node->mod = node->x + node->get_leftmost_child()->x;
    //     }
    // }
    // else {
    //     int left_child_x = node->get_leftmost_child()->x;
    //     int right_child_x = node->get_rightmost_child()->x;
    //     int mid = (left_child_x + right_child_x) / 2;
    //
    //     if(node->is_leftmost()) {
    //         node->x = mid;
    //     }
    //     else {
    //         node->x = node->get_left_sibling()->x + 1;
    //         node->mod = node->x - mid;
    //     }
    // }
    // node->y = depth;
    //
    // if (!node->children.empty() && !node->is_leftmost()) {
    //     check_for_conflicts(node);
    // }
// }

