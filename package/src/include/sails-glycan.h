//
// Created by Jordan Dialpuri on 06/07/2024.
//

#ifndef SAILS_SAILS_GLYCAN_H
#define SAILS_SAILS_GLYCAN_H

#include <list>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <queue>
#include <fstream>

#include "gemmi/model.hpp"

#include "sails-model.h"

namespace Sails {

    struct Sugar {

        Sugar(const std::string &atom, int seqId, Sails::Glycosite& site ) : atom(atom), seqId(seqId), site(site) {}

        std::string atom;
        int seqId;
        Sails::Glycosite site;

        bool operator< (const Sugar& rhs) const
        {
            return seqId < rhs.seqId;
        }
    };

    struct Glycan {
        void add_linkage(int sugar_1, int sugar_2) {
            adjacency_list[sugars[sugar_1].get()].insert(sugars[sugar_2].get());
        }

        void add_sugar(const std::string &atom, int seqId, Sails::Glycosite& residue) {

            if (sugars.find(seqId) != sugars.end()) { // this sugar was already added, don't overwrite
                return;
            }

            sugars[seqId] = std::make_unique<Sugar>(atom, seqId, residue);
        }

        void print_list(gemmi::Structure& structure) {
            for (const auto& pair: adjacency_list) {
                const Sugar* node = pair.first;
                const std::set<Sugar*>& siblings = pair.second;

                gemmi::Residue r = structure.models[node->site.model_idx].chains[node->site.chain_idx].residues[node->site.residue_idx];
                std::cout << r.name << "-" << r.seqid.str() << "-" << node->atom << ": ";

                for (const Sugar* sibling: siblings) {
                    gemmi::Residue r2 = structure.models[sibling->site.model_idx].chains[sibling->site.chain_idx].residues[sibling->site.residue_idx];
                    std::cout << r2.name << "-" << r2.seqid.str()  << "-" << sibling->atom << ", " ;
                }
                std::cout << std::endl;
            }
        }

        std::string get_dot_string(gemmi::Structure& structure, Sails::ResidueDatabase& database) {
            std::string dot = "graph Glycan {\n";
            dot += "rankdir=RL;\n";
            dot += "{\n";

            for (const auto& [root, sugar]: sugars) {
                gemmi::Residue r2 = structure.models[sugar->site.model_idx].chains[sugar->site.chain_idx].residues[sugar->site.residue_idx];
                dot += std::to_string(root);
                dot += " [fillcolor=\"" ;
                dot += database[r2.name].snfg_colour;
                dot += "\" shape=";
                dot += database[r2.name].snfg_shape;
                dot += " style=filled label=\"\"]\n";
            }
            dot += "}\n";

            for (const auto& [root, children]: adjacency_list) {
                for (const auto& child: children) {
                    dot += std::to_string(root->seqId);
                    dot += " -- ";
                    dot += std::to_string(child->seqId);
                    dot += ";\n";
                }
            }

            dot += "}";
            return dot;
        }

        void write_dot_file(const std::string& path, gemmi::Structure& structure, Sails::ResidueDatabase& database) {
            std::ofstream of(path);
            of << get_dot_string(structure, database) << std::endl;
            of.close();
        }

        void bfs(Sugar* root) {
            std::queue<Sugar*> to_visit({root});
            std::set<Sugar*> visited = {root};
            std::map<Sugar*, int> level;
            level[root] = 0;

            while(!to_visit.empty()) {
                Sugar* current_sugar = to_visit.front();
                to_visit.pop();
                int current_level = level[current_sugar];

                for (Sugar* sibling: adjacency_list[current_sugar]) {
                    if (visited.find(sibling) == visited.end()) {  // if we haven't visited this sibling yet
                        to_visit.push(sibling);
                        visited.insert(sibling);
                        level[sibling] = current_level + 1;
                    }
                }
            }
        }

        std::vector<Sugar*> get_terminal_sugars(Sugar* root) {
            std::vector<Sugar*> terminal_sugars;
            dfs(root, terminal_sugars);
            return terminal_sugars;
        }

        void dfs(Sugar* current_sugar, std::vector<Sugar*>& terminal_sugars) {
            std::set<Sugar*>& sugars = adjacency_list[current_sugar];

            if (sugars.empty()) {
                terminal_sugars.push_back( current_sugar );
            }

            for (Sugar* sugar: sugars) {
                dfs(sugar, terminal_sugars);
            }
        }

    private:
        std::map<Sugar*, std::set<Sugar*>> adjacency_list;
        std::map<int, std::unique_ptr<Sugar>> sugars; // used to store sugars until Glycan goes out of scope
    };



}

#endif //SAILS_SAILS_GLYCAN_H
