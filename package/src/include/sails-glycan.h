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
#include <optional>

#include <gemmi/model.hpp>

#include "sails-model.h"
#include "sails-utils.h"

namespace Sails {
    struct Sugar; // forward declaration

    /**
     * @struct Linkage
     * @brief Represents a linkage between two sugar objects.
     *
     * The Linkage struct stores information about a linkage between a donor sugar and an acceptor sugar. It includes
     * pointers to the donor and acceptor sugar, as well as the type of atom in each sugar that is involved in the linkage.
     *
     * @param donor_sugar Pointer to the donor sugar molecule.
     * @param acceptor_sugar Pointer to the acceptor sugar molecule.
     * @param donor_atom The type of atom in the donor sugar molecule that is involved in the linkage.
     * @param acceptor_atom The type of atom in the acceptor sugar molecule that is involved in the linkage.
     */
    struct Linkage {
        Linkage(Sugar *donor_sugar, Sugar *acceptor_sugar, const std::string &donor_atom,
                const std::string &acceptor_atom)
            : donor_sugar(donor_sugar),
              acceptor_sugar(acceptor_sugar),
              donor_atom(donor_atom),
              acceptor_atom(acceptor_atom) {
        }

        Sugar *donor_sugar;
        Sugar *acceptor_sugar;
        std::string donor_atom;
        std::string acceptor_atom;
    };


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
        std::vector<Linkage> linkages;

        bool operator<(const Sugar &rhs) const {
            return seqId < rhs.seqId;
        }

        bool operator==(const Sugar &other) const {
            return this->site == other.site;
        }

        Linkage* find_linkage(const Sugar *sugar1, const Sugar *sugar2) {
            for (auto &linkage: linkages) {
                if (sugar1 == linkage.donor_sugar && sugar2 == linkage.acceptor_sugar) {
                    return &linkage;
                }
            }
            return nullptr;
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

        Glycan(const Glycan &other) : m_structure(other.m_structure),
                                      m_database(other.m_database),
                                      glycosite(other.glycosite) {
            for (const auto &[key, sugar]: other.sugars) {
                sugars[key] = std::make_unique<Sugar>(*sugar);
            }

            for (const auto &[sugar, linkedSugars]: other.adjacency_list) {
                for (auto &linkedSugar: linkedSugars) {
                    adjacency_list[sugars[sugar->site].get()].insert(sugars[linkedSugar->site].get());
                }
            }

            for (const auto &linkage: other.linkage_list) {
                linkage_list.emplace_back(linkage);
            }
        }

        Glycan(gemmi::Structure *structure, ResidueDatabase &database, Glycosite &glycosite): m_structure(structure),
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
        [[nodiscard]] std::map<Glycosite, std::unique_ptr<Sugar> >::const_iterator begin() const {
            return sugars.begin();
        }

        /**
         * @brief Returns a iterator to the end of the map.
         *
         * Returns a constant iterator that points to the last element of the map. This iterator can be used to
         * compare with other iterators to check if it has reached the end of the map.
         *
         * @return A constant iterator pointing to the end of the map.
         * @note The map must not be modified while iterating over it using this iterator.
         */
        [[nodiscard]] std::map<Glycosite, std::unique_ptr<Sugar> >::const_iterator end() const { return sugars.end(); }


        /**
         * @brief Checks if the container is empty.
         *
         * The `empty` method checks if the container is empty. It returns `true` if the container doesn't contain any elements,
         * and `false` otherwise.
         *
         * @return `true` if the container is empty, `false` otherwise.
         */
        [[nodiscard]] bool empty() const {
            return sugars.size() == 1;
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
         * @param donor_atom
         * @param acceptor_atom
         *
         * @note The specified sugar keys must be valid IDs of existing sugar objects.
         *
         * @return void

        */
        void add_linkage(Glycosite &sugar_1, Glycosite &sugar_2, const std::string &donor_atom,
                         const std::string &acceptor_atom) {
            if (sugars.find(sugar_1) == sugars.end()) {
                throw std::runtime_error("Attempted to link a sugar not in the glycan");
            }
            if (sugars.find(sugar_2) == sugars.end()) {
                throw std::runtime_error("Attempted to link a sugar not in the glycan");
            }

            adjacency_list[sugars[sugar_1].get()].insert(sugars[sugar_2].get());

            Linkage linkage = {
                sugars[sugar_1].get(),
                sugars[sugar_2].get(),
                donor_atom,
                acceptor_atom
            };
            linkage_list.emplace_back(linkage);

            // add linkage ptr to the sugars, so we can search the linkages and get information about donor/acceptor atoms
            // make a copy on purpose
            sugars[sugar_1]->linkages.emplace_back(linkage);
            sugars[sugar_2]->linkages.emplace_back(linkage);
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
            if (sugars.find(residue) != sugars.end()) {
                // this sugar was already added, don't overwrite
                return;
            }

            sugars[residue] = std::make_unique<Sugar>(atom, seqId, residue);
        }

        /**
         * @fn Sugar* remove_sugar(int seqId, Sugar* sugar)
         * @brief Remove a sugar molecule from the structure and update the adjacency list.
         *
         * The remove_sugar function removes a sugar molecule with the specified sequence ID from the structure and updates the adjacency list accordingly.
         * It returns a pointer to the sugar molecule that was linked to the removed sugar molecule.
         *
         * @param sugar A pointer to the sugar molecule to be removed.
         * @param erase_from_structure Whether to erase the corresponding residue from the structure member.
         * @return A pointer to the sugar molecule that was linked to the removed sugar molecule, or nullptr if no sugar molecule was linked.
         */
        Sugar *remove_sugar(Sugar *sugar, bool erase_from_structure = true) {
            // erase residue from structure member
            if (erase_from_structure) {
                const auto residue_ptr = &m_structure->models[sugar->site.model_idx].chains[sugar->site.chain_idx].
                        residues;
                residue_ptr->erase(residue_ptr->begin() + sugar->site.residue_idx);
            }

            Sugar *linked_donor = nullptr;
            // remove from adjacency list
            for (auto &[donor, acceptor]: adjacency_list) {
                if (acceptor.find(sugar) != acceptor.end()) {
                    linked_donor = donor;
                    acceptor.erase(sugar);
                }
            }

            // erase seqid from map -> Sugar* will be released so must be done last
            if (sugars.find(sugar->site) != sugars.end()) sugars.erase(sugar->site);

            // return the sugar that was linked (the new terminal sugar)
            return linked_donor;
        }

        /**
         * @brief Finds the previous sugar molecule linked to the given sugar molecule.
         *
         * The find_previous_sugar method searches for the previous sugar molecule linked to the given sugar molecule within
         * the adjacency list. It iterates through the adjacency list and checks if the given sugar molecule is found as an
         * acceptor in any of the pair values. If found, it returns the corresponding donor sugar molecule. If no previous
         * sugar is found, it returns std::nullopt.
         *
         * @param sugar A pointer to the sugar molecule for which to find the previous sugar molecule.
         * @return An optional pointer to the previous sugar molecule linked to the given sugar, or std::nullopt if not found.
         */
        std::optional<Sugar *> find_previous_sugar(Sugar *sugar) const {
            for (auto &[donor, acceptor]: adjacency_list) {
                if (acceptor.find(sugar) != acceptor.end()) {
                    return donor;
                }
            }
            return std::nullopt;
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
        void print_list() const;

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
        [[maybe_unused]] void bfs(Sugar *root);

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
        std::vector<Sugar *> get_terminal_sugars(Glycosite &root_seq_id);

        /**
         * Performs a depth-first search (DFS) on a graph of sugar molecules, starting from
         * a given sugar and collecting terminal sugars.
         *
         * @param current_sugar - The current sugar molecule being visited.
         * @param terminal_sugars - A vector to store the terminal sugar molecules found.
         * @param depth - The depth of the current search
         */
        [[maybe_unused]] void dfs(Sugar *current_sugar, std::vector<Sugar *> &terminal_sugars, int depth);


        /**
         * @brief Get the structure associated with the glycan.
         *
         * This method returns the gemmi::Structure object associated with the glycan.
         *
         * @return The gemmi::Structure object associated with the glycan.
         *
         * @see gemmi::Structure
         */
        [[nodiscard]] gemmi::Structure get_structure() const { return *m_structure; };


        /**
         * @brief Subtract operator for Glycan objects.
         *
         * The operator- subtracts the glycosites of two Glycan objects and returns the set of unique glycosites that exist in the
         * first Glycan object but not in the second Glycan object.
         *
         * @param glycan The Glycan object to subtract from the current Glycan object.
         * @return The set of unique glycosites that exist in the first Glycan object but not in the second Glycan object.
         * @note The function is marked with [[nodiscard]] to indicate that the return value should not be ignored.
         */
        [[nodiscard]] std::set<Glycosite> operator-(const Glycan &glycan);

        /**
         * Internal glycosite [[maybe_unused]]
         */
        Glycosite glycosite;

        // private:

        /**
         * Adjaceny list describing links between sugar pointers
         */
        std::map<Sugar *, std::set<Sugar *> > adjacency_list;

        /**
         * List of linkages used when creating link records
         */
        std::vector<Linkage> linkage_list;


        /**
         * Sugars map used to store sugars in scope until the Glycan is destroyed
         */
        std::map<Glycosite, std::unique_ptr<Sugar> > sugars; // used to store sugars until Glycan goes out of scope

        /**
         * Internal structure pointer - used so that the structure is updated with the caller's structure
         */
        gemmi::Structure *m_structure;

        /**
         * Internal residue database - used to write local dot files [[DEPRECATED]]
         */
        ResidueDatabase m_database;
    };
}

#endif //SAILS_SAILS_GLYCAN_H
