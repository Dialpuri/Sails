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


namespace Sails {

    struct Sugar {

        Sugar(const std::string &atom, int seqId) : atom(atom), seqId(seqId) {}

        std::string atom;
        int seqId;

        bool operator< (const Sugar& rhs) const
        {
            return seqId < rhs.seqId;
        }
    };

    struct Glycan {
        void addLinkage(Sugar& sugar_1, Sugar& sugar_2) {
            unique_nodes.insert(sugar_1);
            unique_nodes.insert(sugar_2);

            adjacency_list[sugar_1].push_back(sugar_2);
        }

        void printGraph() {
            for (const auto& pair : adjacency_list) {
                const Sugar& node = pair.first;
                const std::list<Sugar>& neighbors = pair.second;
                for (const Sugar& neighbor : neighbors) {
                    std::cout << "Root: " << node.seqId << "->" << neighbor.seqId << std::endl;
                }
            }
        }

        void print_list() {
            for (const auto& pair: adjacency_list) {
                const Sugar& node = pair.first;
                const std::list<Sugar>& siblings = pair.second;
                std::cout << node.seqId << ": ";
                for (const Sugar& sibling: siblings) {
                    std::cout << sibling.seqId << ", ";
                }
                std::cout << std::endl;
            }
        }

        void bfs(const Sugar& root) {
            std::queue<Sugar> to_visit({root});
            std::set<Sugar> visited = {root};
            std::map<Sugar, int> level;
            level[root] = 0;

            while(!to_visit.empty()) {
                Sugar current_sugar = to_visit.front();
                to_visit.pop();
                int current_level = level[current_sugar];

                for (const Sugar& sibling: adjacency_list[current_sugar]) {
                    if (visited.find(sibling) == visited.end()) {  // if we haven't visited this sibling yet
                        to_visit.push(sibling);
                        visited.insert(sibling);
                        level[sibling] = current_level + 1;
                    }
                }
            }
        }

        std::vector<Sugar*> get_terminal_sugars(Sugar& root) {
            std::vector<Sugar*> terminal_sugars;
            dfs(root, terminal_sugars);
            return terminal_sugars;
        }

        void dfs(Sugar& current_sugar, std::vector<Sugar*>& terminal_sugars) {
            std::list<Sugar>& sugars = adjacency_list[current_sugar];

            if (sugars.empty()) {
                terminal_sugars.push_back( &current_sugar );
            }

            for (Sugar& sugar: sugars) {
                dfs(sugar, terminal_sugars);
            }
        }

    private:
        std::map<Sugar, std::list<Sugar>> adjacency_list;
        std::set<Sugar> unique_nodes;

    };



}

#endif //SAILS_SAILS_GLYCAN_H
