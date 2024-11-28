//
// Created by jordan on 31/10/24.
//

#include "../include/sails-morph.h"
#include "../include/sails-linkage.h"


void Sails::Morph::transform(Glycan &glycan, PseudoGlycan &pseudo_glycan) {
    std::cout << "Transforming!" << std::endl;

    if (glycan.linkage_count() != pseudo_glycan.linkage_count()) {
        throw std::runtime_error("Requested WURCS must represent a glycan of the same size");
    }

    if (!check_graph_isomorphism(glycan, pseudo_glycan)) {
        throw std::runtime_error("Requested WURCS must have the same structure as the glycan");
    }

    swap_sugars(glycan, pseudo_glycan);
}

void Sails::Morph::swap_sugars(Glycan &glycan, PseudoGlycan &pseudo_glycan) const {
    std::stack<std::pair<Sugar *, Sugar *> > glycan_sugars;
    std::unordered_set<Sugar *> glycan_visited;
    glycan_sugars.emplace(glycan.root_sugar, pseudo_glycan.root_sugar);

    while (!glycan_sugars.empty()) {
        auto [glycan_node, pseudo_glycan_node] = glycan_sugars.top();

        glycan_sugars.pop();

        if (glycan_visited.find(glycan_node) == glycan_visited.end() || glycan_visited.find(pseudo_glycan_node) ==
            glycan_visited.end()) {
            glycan_visited.insert(glycan_node);
            glycan_visited.insert(pseudo_glycan_node);
            }

        std::set<Sugar *> glycan_children = glycan.adjacency_list[glycan_node];
        std::set<Sugar *> pseudo_glycan_children = pseudo_glycan.adjacency_list[pseudo_glycan_node];

        auto glycan_children_iter = glycan_children.begin();
        auto pseudo_glycan_children_iter = pseudo_glycan_children.begin();

        while (glycan_children_iter != glycan_children.end() && pseudo_glycan_children_iter != pseudo_glycan_children.end()) {

            gemmi::Residue* glycan_residue = Utils::get_residue_ptr_from_glycosite((*glycan_children_iter)->site, m_structure);
            gemmi::Residue* pseudo_glycan_residue = Utils::get_residue_ptr_from_glycosite((*pseudo_glycan_children_iter)->site, m_structure);

            if (glycan_residue->name != pseudo_glycan_residue->name) {
                gemmi::Residue replaced_residue = Model::replace_residue(glycan_residue, pseudo_glycan_residue->name);
                gemmi::Chain* chain_ptr = Utils::get_chain_ptr_from_glycosite((*glycan_children_iter)->site, m_structure);
                chain_ptr->residues[(*glycan_children_iter)->site.residue_idx] = replaced_residue;
            }

            glycan_sugars.emplace(*glycan_children_iter, *pseudo_glycan_children_iter);
            ++glycan_children_iter;
            ++pseudo_glycan_children_iter;
        }
    }
}

bool Sails::Morph::check_graph_isomorphism(Glycan &glycan, PseudoGlycan &pseudo_glycan) {
    std::stack<std::pair<Sugar *, Sugar *> > glycan_sugars;

    std::unordered_set<Sugar *> glycan_visited;
    bool isomorphous = true;

    glycan_sugars.emplace(glycan.root_sugar, pseudo_glycan.root_sugar);

    while (!glycan_sugars.empty()) {
        auto [glycan_node, pseudo_glycan_node] = glycan_sugars.top();

        glycan_sugars.pop();

        if (glycan_visited.find(glycan_node) == glycan_visited.end() || glycan_visited.find(pseudo_glycan_node) ==
            glycan_visited.end()) {
            glycan_visited.insert(glycan_node);
            glycan_visited.insert(pseudo_glycan_node);
        }

        std::set<Sugar *> glycan_children = glycan.adjacency_list[glycan_node];
        std::set<Sugar *> pseudo_glycan_children = pseudo_glycan.adjacency_list[pseudo_glycan_node];

        if (glycan_children.size() != pseudo_glycan_children.size()) {
            isomorphous = false;
            break;
        }

        auto glycan_children_iter = glycan_children.begin();
        auto pseudo_glycan_children_iter = pseudo_glycan_children.begin();

        while (glycan_children_iter != glycan_children.end() && pseudo_glycan_children_iter != pseudo_glycan_children.end()) {

            glycan_sugars.emplace(*glycan_children_iter, *pseudo_glycan_children_iter);

            ++glycan_children_iter;
            ++pseudo_glycan_children_iter;
        }
    }

    return isomorphous;
}
