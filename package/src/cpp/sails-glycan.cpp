//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../include/sails-glycan.h"

void Sails::Glycan::print_list() {
    for (const auto &pair: adjacency_list) {
        const Sugar *node = pair.first;
        const std::set<Sugar *> &siblings = pair.second;

        gemmi::Residue r = m_structure.models[node->site.model_idx].chains[node->site.chain_idx].residues[node->site.residue_idx];
        std::cout << r.name << "-" << r.seqid.str() << "-" << node->atom << ": ";

        for (const Sugar *sibling: siblings) {
            gemmi::Residue r2 = m_structure.models[sibling->site.model_idx].chains[sibling->site.chain_idx].residues[sibling->site.residue_idx];
            std::cout << r2.name << "-" << r2.seqid.str() << "-" << sibling->atom << ", ";
        }
        std::cout << std::endl;
    }
}

void Sails::Glycan::print_sugars() {
    for (const auto& [index, sugar]: sugars) {
        std::cout << "Seqid " << index << std::endl;
    }
}

std::string Sails::Glycan::get_dot_string() {
    std::string dot = "graph Glycan {\n";
    dot += "rankdir=RL;\n";
    dot += "{\n";

    for (const auto &[root, sugar]: sugars) {
        gemmi::Residue r2 = m_structure.models[sugar->site.model_idx].chains[sugar->site.chain_idx].residues[sugar->site.residue_idx];
        dot += std::to_string(root);
        dot += " [fillcolor=\"";
        dot += m_database[r2.name].snfg_colour;
        dot += "\" shape=";
        dot += m_database[r2.name].snfg_shape;
        dot += " style=filled label=\"\"]\n";
    }
    dot += "}\n";

    for (const auto &[root, children]: adjacency_list) {
        for (const auto &child: children) {
            dot += std::to_string(root->seqId);
            dot += " -- ";
            dot += std::to_string(child->seqId);
            dot += ";\n";
        }
    }

    dot += "}";
    return dot;
}

void
Sails::Glycan::write_dot_file(const std::string &path) {
    std::ofstream of(path);
    of << get_dot_string() << std::endl;
    of.close();
}

void Sails::Glycan::bfs(Sails::Sugar *root) {
    std::queue<Sugar *> to_visit({root});
    std::set<Sugar *> visited = {root};
    std::map<Sugar *, int> level;
    level[root] = 0;

    while (!to_visit.empty()) {
        Sugar *current_sugar = to_visit.front();
        to_visit.pop();
        int current_level = level[current_sugar];

        for (Sugar *sibling: adjacency_list[current_sugar]) {
            if (visited.find(sibling) == visited.end()) {  // if we haven't visited this sibling yet
                to_visit.push(sibling);
                visited.insert(sibling);
                level[sibling] = current_level + 1;
            }
        }
    }
}

void Sails::Glycan::dfs(Sails::Sugar *current_sugar, std::vector<Sugar *> &terminal_sugars) {
    std::set<Sugar *> &sugar_set = adjacency_list[current_sugar];

    if (sugar_set.empty()) {
        terminal_sugars.push_back(current_sugar);
    }

    for (Sugar *sugar: sugar_set) {
        dfs(sugar, terminal_sugars);
    }
}

std::vector<Sails::Sugar *> Sails::Glycan::get_terminal_sugars(int root_seq_id) {
    {
        if (sugars.find(root_seq_id) == sugars.end()) {
            throw std::runtime_error("Root SeqId is not valid : " + std::to_string(root_seq_id));
        }
        std::vector<Sugar *> terminal_sugars;
        dfs(sugars[root_seq_id].get(), terminal_sugars);
        return terminal_sugars;
    }
}