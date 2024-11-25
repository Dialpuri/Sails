//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../include/sails-glycan.h"


std::vector<Sails::Glycosite> Sails::Glycan::find_children(Sugar *sugar) {
    std::stack<Sugar *> sugar_stack;
    std::set<Sugar *> visited;
    sugar_stack.push(sugar);

    std::vector<Glycosite> child_sites;
    while (!sugar_stack.empty()) {
        Sugar *current_sugar = sugar_stack.top();
        sugar_stack.pop();

        if (visited.find(current_sugar) == visited.end()) {
            if (current_sugar != sugar) child_sites.emplace_back(current_sugar->site);
            visited.insert(current_sugar);
        }

        for (auto &child: adjacency_list[current_sugar]) {
            if (visited.find(child) == visited.end()) {
                sugar_stack.push(child);
            }
        }
    }
    return child_sites;
}

std::vector<Sails::LinkageData> Sails::Glycan::extract_linkage_data(LinkageDatabase &linkage_database) const {
    std::vector<LinkageData> linkages;
    for (auto &linkage: linkage_list) {
        const gemmi::Residue *donor_residue = Sails::Utils::get_residue_ptr_from_glycosite(
            linkage.donor_sugar->site, m_structure);
        gemmi::Residue *acceptor_residue = Sails::Utils::get_residue_ptr_from_glycosite(
            linkage.acceptor_sugar->site, m_structure);

        std::vector<LinkageData> available_linkages = linkage_database[donor_residue->name];

        auto found_linkage = std::find_if(available_linkages.begin(), available_linkages.end(),
                                          [&](const Sails::LinkageData &available_linkage) {
                                              return available_linkage.donor_number == linkage.donor_number &&
                                                     available_linkage.acceptor == acceptor_residue->name;
                                          });
        if (found_linkage != available_linkages.end()) {
            linkages.emplace_back(*found_linkage);
        }
    }

    if (linkages.size() != linkage_list.size()) {
        throw std::runtime_error("Linkage list and extracted linkages are not the same size, check data file.");
    }
    return linkages;
}

void Sails::Glycan::print_list() const {
    std::cout << "Adjacency List" << std::endl;
    for (const auto &pair: adjacency_list) {
        const Sugar *node = pair.first;
        const std::set<Sugar *> &siblings = pair.second;

        gemmi::Residue r = m_structure->models[node->site.model_idx].chains[node->site.chain_idx].residues[node->site.
            residue_idx];
        std::cout << r.name << "-" << r.seqid.str() << "-" << node->atom << ": ";

        for (const Sugar *sibling: siblings) {
            gemmi::Residue r2 = m_structure->models[sibling->site.model_idx].chains[sibling->site.chain_idx].residues[
                sibling->site.residue_idx];
            std::cout << r2.name << "-" << r2.seqid.str() << "-" << sibling->atom << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Sails::Glycan::print_sugars() {
    for (const auto &[index, sugar]: sugars) {
        std::cout << "Sites " << Utils::format_site_key(index) << std::endl;
    }
}

std::string Sails::Glycan::get_dot_string() {
    std::string dot = "graph Glycan {\n";
    dot += "rankdir=RL;\n";
    dot += "{\n";

    for (const auto &[root, sugar]: sugars) {
        gemmi::Residue r2 = m_structure->models[sugar->site.model_idx].chains[sugar->site.chain_idx].residues[sugar->
            site.residue_idx];
        dot += Utils::format_site_key(sugar->site);
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

void Sails::Glycan::write_dot_file(const std::string &path) {
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
            if (visited.find(sibling) == visited.end()) {
                // if we haven't visited this sibling yet
                to_visit.push(sibling);
                visited.insert(sibling);
                level[sibling] = current_level + 1;
            }
        }
    }
}

void Sails::Glycan::dfs(Sugar *current_sugar, std::vector<Sugar *> &terminal_sugars, int depth = 0) {
    std::set<Sugar *> &sugar_set = adjacency_list[current_sugar];
    if (sugar_set.empty()) {
        current_sugar->depth = depth;
        terminal_sugars.push_back(current_sugar);
    }

    for (Sugar *sugar: sugar_set) {
        sugar->depth = depth + 1;
        dfs(sugar, terminal_sugars, depth + 1);
    }
}

void Sails::Glycan::dfs_sites(Sugar *current_sugar, std::vector<Glycosite> &sites, int depth) {
    const std::set<Sugar *> &sugar_set = adjacency_list[current_sugar];
    current_sugar->depth = depth;
    sites.push_back(current_sugar->site);

    for (Sugar *sugar: sugar_set) {
        sugar->depth = depth + 1;
        dfs_sites(sugar, sites, depth + 1);
    }
}

std::set<Sails::Glycosite> Sails::Glycan::operator-(const Glycan &glycan) {
    std::set<Glycosite> this_keys;
    std::transform(this->sugars.begin(), this->sugars.end(), std::inserter(this_keys, this_keys.end()),
                   [](auto &kv_pair) { return kv_pair.first; });

    std::set<Glycosite> other_keys;
    std::transform(glycan.sugars.begin(), glycan.sugars.end(), std::inserter(other_keys, other_keys.end()),
                   [](auto &kv_pair) { return kv_pair.first; });

    std::set<Glycosite> difference;
    std::set_difference(this_keys.begin(), this_keys.end(), other_keys.begin(), other_keys.end(),
                        std::inserter(difference, difference.begin()));
    return difference;
}

std::vector<Sails::Sugar *> Sails::Glycan::get_terminal_sugars(Glycosite &root_seq_id) {
    if (sugars.find(root_seq_id) == sugars.end()) {
        throw std::runtime_error("Root SeqId is not valid");
    }
    std::vector<Sugar *> terminal_sugars;
    dfs(sugars[root_seq_id].get(), terminal_sugars);
    return terminal_sugars;
}
