//
// Created by Jordan Dialpuri on 04/08/2024.
//

#include "../../include/snfg/sails-snfg.h"

#include "src/include/snfg/sails-snfg-shape.h"

std::string Sails::SNFG::create_svg_header() const {
    return "<svg width=\"" + std::to_string(SVG_WIDTH) + "\" height=\"" + std::to_string(SVG_HEIGHT) +
           "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
}

std::string Sails::SNFG::create_svg_footer() {
    return "</svg>\n";
}

Sails::SVGStringObject Sails::SNFG::create_svg_line(int x1, int y1, int x2, int y2) {
    std::string object =
            "<line x1=\"" + std::to_string(x1) + "\" y1=\"" + std::to_string(y1) + "\" x2=\"" + std::to_string(x2) +
            "\" y2=\"" + std::to_string(y2) + "\" stroke=\"black\" stroke-width=\"4\" />\n";
    return {object, 2};
}

Sails::SVGStringObject Sails::SNFG::create_svg_text(int x, int y, const std::string &text) {
    std::string object = "<text x=\"" + std::to_string(x) + "\" y=\"" + std::to_string(y) +
                         "\" font-family=\"Verdana\" font-size=\"20\" fill=\"black\" text-anchor=\"middle\">" + text +
                         "</text>\n";
    return {object, 0};
}

Sails::SVGStringObject Sails::SNFG::create_donor_labels(SNFGNode *parent, SNFGNode *node, Linkage *linkage) {

    int dx = node->x - parent->x;
    int dy = node->y - parent->y;

    double ratio = 0.3;
    double xt = parent->x;
    double yt = parent->y;
    if (dx != 0) {
        int xbuffer = -10;
        int ybuffer = dy == 0 ? 20: 19;
        ybuffer += dy > 0 ? (dy*0.2) : 0;
        xt += (ratio*dx) + xbuffer;
        yt += (ratio*dy) + ybuffer;
    }
    else {
        yt += (ratio*dy) + 20;
        xt += 10;
    }
    std::string donor_number = {linkage->donor_atom[linkage->donor_atom.size()-1]};
    return create_svg_text(xt, yt, donor_number);
}

Sails::SVGStringObject Sails::SNFG::create_anomer_labels(SNFGNode* parent, SNFGNode* node, Linkage* linkage) {
    int dx = node->x - parent->x;
    int dy = node->y - parent->y;

    double ratio = 0.7;

    double xt = parent->x;
    double yt = parent->y;
    if (dx != 0) {
        int xbuffer = dy > 0 ? 10: 5;
        int ybuffer = dy == 0 ? 20: 19;
        ybuffer += dy > 0 ? (dy*0.2) : 0;
        xt += (ratio*dx) + xbuffer;
        yt += (ratio*dy) + ybuffer;
    }
    else {
        yt += (ratio*dy) + 10;
        xt += 10;
    }


    gemmi::Residue* acceptor_residue_ptr = Utils::get_residue_ptr_from_glycosite(node->sugar->site, m_structure);
    std::string acceptor_name = acceptor_residue_ptr->name;
    std::string donor_number = {m_database->operator[](acceptor_name).anomer};
    return create_svg_text(xt, yt, donor_number);
}

void Sails::SNFG::printTree(SNFGNode *root, Sails::SNFGNode *node, int level) {
    if (node == nullptr) return;

    std::cout << "Level: " << level << " -> node: "
            << Sails::Utils::format_residue_from_site(node->sugar->site, m_structure) << ", x: " << node->y << ", y: "
            << node->x << std::endl;

    for (auto &child: node->children) {
        printTree(root, child.get(), level + 1);
    }
}


void Sails::SNFG::check_for_conflicts(SNFGNode *node) {
    int min_distance = tree_distance + node_size;
    int shift = 0;

    std::map<int, int> node_contour;
    get_lcontour(node, 0, node_contour);

    auto sibling = node->get_leftmost_sibling();
    while (sibling != nullptr && sibling != node) {
        std::map<int, int> sibling_contour;
        get_rcontour(sibling, 0, sibling_contour);

        auto last_sibling_contour = std::prev(sibling_contour.end())->first;
        auto last_node_contour = std::prev(node_contour.end())->first;
        for (int i = node->x + 1; i <= std::min(last_sibling_contour, last_node_contour); i++) {
            int distance = node_contour[i] - sibling_contour[i];
            if (distance + shift < min_distance) {
                shift = min_distance - distance;
            }
        }

        if (shift > 0) {
            node->y += shift;
            node->mod += shift;

            center_nodes_between(node, sibling);
            shift = 0;
        }
        sibling = sibling->get_right_sibling();
    }
}

void Sails::SNFG::center_nodes_between(SNFGNode *lnode, SNFGNode *rnode) {
    int li = lnode->parent_position;
    int ri = rnode->parent_position;

    int node_between = (ri - li) - 1;
    if (node_between > 0) {
        int node_distance = (lnode->y - rnode->y) / (node_between + 1);
        int count = 1;
        for (int i = li + 1; i < ri; i++) {
            SNFGNode *mid_node = lnode->parent->children[i].get();
            int desired_x = rnode->y + (node_distance * count);
            int offset = desired_x - mid_node->y;
            mid_node->y += offset;
            mid_node->y += offset;
            count++;
        }
        check_for_conflicts(lnode);
    }
}


void Sails::SNFG::get_lcontour(Sails::SNFGNode *node, int mod_sum, std::map<int, int> &values) {
    if (values.find(node->x) == values.end()) {
        values[node->x] = node->y + mod_sum;
    } else {
        values[node->x] = std::min(values[node->x], node->y + mod_sum);
    }

    mod_sum += node->mod;

    for (auto &child: node->children) {
        get_lcontour(child.get(), mod_sum, values);
    }
}

void Sails::SNFG::get_rcontour(Sails::SNFGNode *node, int mod_sum, std::map<int, int> &values) {
    if (values.find(node->x) == values.end()) {
        values[node->x] = node->y + mod_sum;
    } else {
        values[node->x] = std::max(values[node->x], node->y + mod_sum);
    }

    mod_sum += node->mod;

    for (auto &child: node->children) {
        get_rcontour(child.get(), mod_sum, values);
    }
}


void Sails::SNFG::calculate_node_positions(Sails::SNFGNode *node) {
    initialise_nodes(node, 0);

    calculate_initial_positions(node);

    check_all_children_on_screen(node);

    calculate_final_positions(node, 0);
}

void Sails::SNFG::initialise_nodes(Sails::SNFGNode *node, int depth) {
    node->y = -1;
    node->x = depth;
    node->mod = 0;

    for (auto &child: node->children) {
        initialise_nodes(child.get(), depth + 1);
    }
}

void Sails::SNFG::calculate_initial_positions(Sails::SNFGNode *node) {
    for (auto &child: node->children) {
        calculate_initial_positions(child.get());
    }

    if (node->is_leaf()) {
        if (!node->is_leftmost()) {
            node->y = node->get_left_sibling()->y + node_size + sibling_distance;
        } else {
            node->y = 0;
        }
    } else if (node->children.size() == 1) {
        if (node->is_leftmost()) {
            node->y = node->get_leftmost_child()->y;
        } else {
            node->y = node->get_left_sibling()->y + node_size + sibling_distance;
            node->mod = node->y - node->get_leftmost_child()->y;
        }
    } else {
        int lc = node->get_leftmost_child()->y;
        int rc = node->get_rightmost_child()->y;

        int mid = (lc + rc) / 2;
        if (node->is_leftmost()) {
            node->y = mid;
        } else {
            node->y = node->get_left_sibling()->y + node_size + sibling_distance;
            node->mod = node->y - mid;
        }
    }
    if (node->children.size() > 0 && !node->is_leftmost()) {
        check_for_conflicts(node);
    }
}


void Sails::SNFG::check_all_children_on_screen(SNFGNode *node) {
    std::map<int, int> node_contour;

    get_lcontour(node, 0, node_contour);

    int shift_amount = 0;
    for (auto &[k, v]: node_contour) {
        if (node_contour[k] + shift_amount < 0) {
            shift_amount = -node_contour[k];
        }
    }
    if (shift_amount > 0) {
        node->y += shift_amount;
        node->mod += shift_amount;
    }
}

void Sails::SNFG::calculate_final_positions(Sails::SNFGNode *node, float mod_sum) {
    node->y += mod_sum;
    mod_sum += node->mod;

    for (auto &child: node->children) {
        calculate_final_positions(child.get(), mod_sum);
    }

    if (node->children.empty()) {
    }

    node->x = SVG_WIDTH - 200 - (100 * node->x);
    node->y += 400;
}


void Sails::SNFG::assign_snfg_types(Sugar *child, SNFGNode *new_child) {
    const gemmi::Residue *residue_ptr = Utils::get_residue_ptr_from_glycosite(child->site, m_structure);
    if (m_database->find(residue_ptr->name) != m_database->end()) {
        const ResidueData *instance = &m_database->operator[](residue_ptr->name);
        new_child->snfg_colour = instance->snfg_colour;
        new_child->snfg_shape = instance->snfg_shape;
        new_child->special = instance->special;
    }
}

void Sails::SNFG::form_snfg_node_system(SNFGNode *root, Sugar *sugar, Glycan &glycan, std::vector<std::unique_ptr<SNFGNode>> &special_nodes) {
    int pos = 0;

    root->residue = Utils::get_residue_ptr_from_glycosite(root->sugar->site, m_structure);
    root->chain = Utils::get_chain_ptr_from_glycosite(root->sugar->site, m_structure);

    // transform set of linked sugars, to a vector, so we can sort by donor_number. This affects what node is positioned
    // where in the SNFG. We want lower donor numbers to be on top, so sort here.
    std::vector<Sugar *> adjacency_list = std::vector(glycan.adjacency_list[sugar].begin(),
                                                      glycan.adjacency_list[sugar].end());
    std::sort(adjacency_list.begin(), adjacency_list.end(), [&](const Sugar* s1, const Sugar* s2) {
        const Linkage *l1 = sugar->find_linkage(sugar, s1);
        const Linkage *l2 = sugar->find_linkage(sugar, s2);
        if (l1 == nullptr || l2 == nullptr) return false;
        return l1->donor_number < l2->donor_number;
    });

    for (auto &child: adjacency_list) {
        auto new_child = std::make_unique<SNFGNode>(child);
        new_child->parent_position = pos;
        new_child->parent = root;

        assign_snfg_types(child, new_child.get());

        if (!new_child->special) {
            root->children.push_back(std::move(new_child));
        } else {
            new_child->residue = Utils::get_residue_ptr_from_glycosite(new_child->sugar->site, m_structure);
            new_child->chain = Utils::get_chain_ptr_from_glycosite(new_child->sugar->site, m_structure);
            special_nodes.push_back(std::move(new_child));
            continue;
        }
        form_snfg_node_system(root->children.back().get(), child, glycan, special_nodes);
        pos++;
    }
}


void Sails::SNFG::create_svg(std::vector<SVGStringObject> &snfg_objects, SNFGNode *parent, SNFGNode *node) {
    // find linkage
    Linkage *linkage = parent->sugar->find_linkage(parent->sugar, node->sugar);
    if (linkage != nullptr) {
        node->linkage = linkage;
        snfg_objects.emplace_back(create_donor_labels(parent, node, linkage));
        snfg_objects.emplace_back(create_anomer_labels(parent, node, linkage));

    }

    // get shape
    snfg_objects.emplace_back(get_svg_shape(node)->draw());

    // get line
    snfg_objects.emplace_back(create_svg_line(parent->x, parent->y, node->x, node->y));
    // snfg_objects.emplace_back(create_svg_text(parent->x, parent->y, linkage->donor_atom));
    for (auto &child: node->children) {
        create_svg(snfg_objects, node, child.get());
    }
}


void Sails::SNFG::order_svg(std::vector<Sails::SVGStringObject> &objects) {
    std::sort(objects.begin(), objects.end(),
              [&](Sails::SVGStringObject &obj1, Sails::SVGStringObject &obj2) {
                  return obj1.priority > obj2.priority;
              });
}

void Sails::SNFG::position_special_nodes(std::vector<std::unique_ptr<SNFGNode>>& nodes) {
    if (nodes.empty()) return;
    for (auto &child: nodes) {
        child->x = child->parent->x;
        child->y = child->parent->y + node_size;
        child->parent->children.emplace_back(std::move(child));
    }
}


std::string Sails::SNFG::create_snfg(Glycan &glycan, Glycosite &base_residue) {
    auto base_sugar = glycan.sugars[base_residue].get();

    auto root = std::make_unique<SNFGNode>(base_sugar);
    std::vector<std::unique_ptr<SNFGNode>> special_nodes;
    form_snfg_node_system(root.get(), base_sugar, glycan, special_nodes);
    calculate_node_positions(root.get());
    position_special_nodes(special_nodes);

    std::vector<Sails::SVGStringObject> objects;
    create_svg(objects, root.get(), root.get());
    order_svg(objects);

    std::stringstream f;
    f << create_svg_header();
    for (const auto &obj: objects) {
        f << obj.object;
    }
    f << create_svg_footer();
    return f.str();
}
