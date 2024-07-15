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
    /**
     * @struct Sugar
     * @brief Represents a sugar molecule.
     *
     * The Sugar struct stores information about a sugar molecule. It includes the type of atom, the sequence ID,
     * and the glycosite of the sugar molecule.
     *
     * @param atom The type of atom in the sugar molecule.
     * @param seqId The sequence ID of the sugar molecule.
     * @param site The glycosite of the sugar molecule.
     */
    struct Sugar {
        Sugar(const Sugar &sugar) {
            atom = sugar.atom;
            seqId = sugar.seqId;
            site = sugar.site;
            depth = sugar.depth;
        }

        Sugar(const std::string &atom, int seqId, Glycosite &site) : atom(atom), seqId(seqId), site(site) {
            depth = 0;
        }

        Sugar(const std::string &atom, int seqId, Glycosite &site, int depth) : atom(atom), seqId(seqId), site(site),
            depth(depth) {
        }

        std::string atom;
        int seqId{};
        Glycosite site;
        int depth{};

        bool operator<(const Sugar &rhs) const {
            return seqId < rhs.seqId;
        }
    };

    /**
     * @brief Glycan represents a glycan structure.
     *
     * The Glycan class is responsible for constructing and managing a glycan structure composed of sugars.
     * It provides methods for adding linkages between sugars, adding sugars to the glycan, printing the glycan
     * structure, generating a DOT string representation of the glycan, writing the glycan to a DOT file,
     * performing breadth-first search (BFS) traversal on the glycan structure, retrieving terminal sugars,
     * and performing depth-first search (DFS) traversal on the glycan structure.
     */
    struct Glycan {
        Glycan() = default;

        Glycan(gemmi::Structure &structure, ResidueDatabase &database, Glycosite &glycosite): m_structure(structure),
            m_database(database), glycosite(glycosite) {
        }


        /**
         * @brief Returns an iterator to the beginning of the map.
         *
         * Returns a constant iterator that points to the first element of the map. This iterator can be used to traverse the map
         * in a forward direction.
         *
         * @return A constant iterator that points to the first element of the map.
         * @note The map must not be modified while iterating over it using this iterator.
         */
        [[nodiscard]] std::map<int, std::unique_ptr<Sugar> >::const_iterator begin() const { return sugars.begin(); }

        /**
         * @brief Returns a iterator to the end of the map.
         *
         * Returns a constant iterator that points to the last element of the map. This iterator can be used to
         * compare with other iterators to check if it has reached the end of the map.
         *
         * @return A constant iterator pointing to the end of the map.
         * @note The map must not be modified while iterating over it using this iterator.
         */
        [[nodiscard]] std::map<int, std::unique_ptr<Sugar> >::const_iterator end() const { return sugars.end(); }


        /**
         * @brief Checks if the container is empty.
         *
         * The `empty` method checks if the container is empty. It returns `true` if the container doesn't contain any elements,
         * and `false` otherwise.
         *
         * @return `true` if the container is empty, `false` otherwise.
         */
        [[nodiscard]] bool empty() const {
            return sugars.empty();
        }


        /**
         * @brief Returns the size of the sugar list.
         *
         * The size of the sugar list represents the number of elements stored in it.
         *
         * @return The size of the sugar list as a `size_t` value.
         */
        [[nodiscard]] size_t size() const {
            return sugars.size();
        }

        /**
         * @brief Adds linkage between two sugars.
         *
         * This function creates a linkage between two sugars by the seqId. Sugars must have been added with add_sugar
         * before this function is called, if not, an error is raised.
         *
         * @param sugar_1_key The ID of the first sugar object.
         * @param sugar_2_key The ID of the second sugar object.
         *
         * @note The specified sugar keys must be valid IDs of existing sugar objects.
         *
         * @return void

        */
        void add_linkage(int sugar_1_key, int sugar_2_key) {
            if (sugars.find(sugar_1_key) == sugars.end()) {
                throw std::runtime_error("Attempted to link a sugar not in the glycan");
            }
            if (sugars.find(sugar_2_key) == sugars.end()) {
                throw std::runtime_error("Attempted to link a sugar not in the glycan");
            }

            adjacency_list[sugars[sugar_1_key].get()].insert(sugars[sugar_2_key].get());
        }

        /**
         * @brief Adds a sugar molecule to be stored in this glycan.
         *
         * The add_sugar method adds a sugar molecule to the collection of sugars. If the sugar molecule with the given seqId
         * already exists in the collection, it will not be overwritten.
         *
         * @param atom The type of atom in the sugar molecule.
         * @param seqId The sequence ID of the sugar molecule.
         * @param residue The glycosite of the sugar molecule.
         */
        void add_sugar(const std::string &atom, int seqId, Glycosite &residue) {
            if (sugars.find(seqId) != sugars.end()) {
                // this sugar was already added, don't overwrite
                return;
            }

            sugars[seqId] = std::make_unique<Sugar>(atom, seqId, residue);
        }

        /**
         * @fn Sugar* remove_sugar(int seqId, Sugar* sugar)
         * @brief Remove a sugar molecule from the structure and update the adjacency list.
         *
         * The remove_sugar function removes a sugar molecule with the specified sequence ID from the structure and updates the adjacency list accordingly.
         * It returns a pointer to the sugar molecule that was linked to the removed sugar molecule.
         *
         * @param sugar A pointer to the sugar molecule to be removed.
         * @return A pointer to the sugar molecule that was linked to the removed sugar molecule, or nullptr if no sugar molecule was linked.
         */
        Sugar* remove_sugar(Sugar* sugar) {

            // erase residue from structure member
            const auto residue_ptr = &m_structure.models[sugar->site.model_idx].chains[sugar->site.chain_idx].residues;
            residue_ptr->erase(residue_ptr->begin()+sugar->site.residue_idx);

            Sugar* linked_donor = nullptr;
            // remove from adjacency list
            for (auto& [donor, acceptor]: adjacency_list) {
                if (acceptor.find(sugar) != acceptor.end()) {
                    linked_donor = donor;
                    acceptor.erase(sugar);
                }
            }

            // erase seqid from map -> Sugar* will be released so must be done last
            sugars.erase(sugar->seqId);

            // return the sugar that was linked (the new terminal sugar)
            return linked_donor;
        }

        /**
        * @brief Prints the adjacency list of the glycan structure.
        *
        * The print_list() function iterates over the adjacency list of the glycan structure and prints each sugar and its linked siblings.
        * It retrieves the residue information from the gemmi::Structure and gemmi::Residue objects and prints them as part of the output.
        * The output is printed to the console.
        *
        * @return void
        */
        void print_list();

        /**
        * @brief Prints the sugars contained within glycan structure.
        *
        * @return void
        */
        void print_sugars();

        /**
         * @brief Retrieves a DOT representation of the glycan structure.
         *
         * The get_dot_string() function generates a DOT string representation of the glycan structure.
         * The function uses the adjacency list to traverse the sugars and their linked siblings to construct the DOT string.
         * The DOT string represents the relationships between sugars as edges and the sugars themselves as nodes.
         *
         * @return A string containing the DOT representation of the glycan structure.
         */
        std::string get_dot_string();

        /**
         * Writes a dot file at the specified path.
         *
         * This function takes a path to a dot file as input and writes a dot file at the specified location.
         * The dot file is used for creating graphs using Graphviz graph visualization software.
         *
         * @param path The path to the dot file to be written.
         */
        void write_dot_file(const std::string &path);

        /**
         * @brief Performs a breadth-first search (BFS) traversal on the glycan structure starting from the specified root sugar.
         *
         * The bfs() function implements the breadth-first search algorithm to traverse the glycan structure in a level-wise manner.
         * It starts from the root sugar and visits all its linked sugars before moving on to the next level of sugars.
         * The function uses a queue to keep track of the sugars to be visited and a set to keep track of visited sugars.
         * It also uses a map to store the level of each sugar in the traversal.
         *
         * @param root A pointer to the root sugar from which the traversal starts.
         *
         * @return void
         */
        void bfs(Sugar *root);

        /**
         * @brief This function retrieves all the terminal sugars in a tree structure starting from the given root sugar.
         *
         * The function takes a pointer to the root sugar of a tree structure as input and returns a collection of all the terminal sugars
         * in that tree structure.
         *
         * @param root The seqId to the root residue of the tree structure.
         * @throws std::runtime_error If supplied root_seq_id is invalid.
         *
         * @return A collection of terminal sugars in the tree structure.
         */
        std::vector<Sugar *> get_terminal_sugars(int root_seq_id);

        /**
         * Performs a depth-first search (DFS) on a graph of sugar molecules, starting from
         * a given sugar and collecting terminal sugars.
         *
         * @param current_sugar - The current sugar molecule being visited.
         * @param terminal_sugars - A vector to store the terminal sugar molecules found.
         * @param depth - The depth of the current search
         */
        void dfs(Sugar *current_sugar, std::vector<Sugar *> &terminal_sugars, int depth);


        /**
         * @brief Get the structure associated with the glycan.
         *
         * This method returns the gemmi::Structure object associated with the glycan.
         *
         * @return The gemmi::Structure object associated with the glycan.
         *
         * @see gemmi::Structure
         */
        [[nodiscard]] gemmi::Structure get_structure() const {return m_structure;};

        Glycosite glycosite;

    private:
        std::map<Sugar *, std::set<Sugar *> > adjacency_list;
        std::map<int, std::unique_ptr<Sugar> > sugars; // used to store sugars until Glycan goes out of scope

        gemmi::Structure m_structure;
        ResidueDatabase m_database;
    };
}

#endif //SAILS_SAILS_GLYCAN_H
