//
// Created by Jordan Dialpuri on 19/02/2024.
//

#include "sails-analyse.h"

#include "sails-data.h"

void printTree(const Node* root, int level = 0) {
    if(root)
    {
        // Print indentation
        std::cout << std::string(level * 4, ' ');

        // Print node data
        std::cout << root->data.get_m() << std::endl;

        // Recursive call for each child
        for(const auto& child : root->children)
        {
            printTree(child.get(), level + 1);
        }
    }
}

bool SailsAnalyse::get_value(clipper::MiniMol& mol, SailsData& database, std::unique_ptr<Node>& node) {
    if (node == nullptr) {return false;}

    clipper::MMonomer residue = node->data.get_mmonomer(mol);
    Sails::ResidueType type = Sails::ResidueMap[residue.type()];
    if (database.data_map.find(type) == database.data_map.end()) { std::cout << "Data on this residue is not available\n"; return true;}

    Data* data = database.data_map[type];

    for (const std::string& donor: data->get_donors()) {
        int atom_index = residue.lookup(donor, clipper::MM::UNIQUE);
        if (atom_index < 0) {std::cout << "Can't find donor atom\n"; continue;}

        clipper::Coord_orth atom_coord_orth = residue[atom_index].coord_orth();
        std::vector<clipper::MAtomIndexSymmetry> atoms = ns.atoms_near(atom_coord_orth, 2.5);

        std::set<NeighbourSearchResult> nearby_residues;
        for (const auto& atom: atoms) {
            nearby_residues.insert({atom.polymer(), atom.monomer()});
        }

        std::set<std::string> possible_links = data->get_possible_links();

        for (const auto& nearby_residue: nearby_residues) {
            if (nearby_residue.is_same(node.get())) {continue;}

            clipper::MMonomer nearby_mmonomer = nearby_residue.extract_monomer(mol);
            if (possible_links.find(nearby_mmonomer.type()) == possible_links.end()) {continue;}

            float min_dist = 1e3;
            int closest_atom_idx = -1;
            for (int a = 0; a < nearby_mmonomer.size(); a++) {
                clipper::Coord_orth co = nearby_mmonomer[a].coord_orth();
                float dist = clipper::Coord_orth::length(atom_coord_orth, co);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest_atom_idx = a;
                };
            }
            if (closest_atom_idx < 0) {std::cout << "Error in closest atom algo."; continue;}

            NodeData node_data = {nearby_residue.p, nearby_residue.m, nearby_mmonomer[closest_atom_idx].name(), donor};
            node->add_child(node_data);
        }
    }

    for (auto& child: node->children) {
        get_value(mol, database, child);
    }

    return false;
}

SailsAnalyse::SailsAnalyse(clipper::MiniMol& mol) {

    ns = clipper::MAtomNonBond(mol, 2);
    std::vector<std::unique_ptr<Node>> nodes;


    SailsData database;

    for (int p = 0; p < mol.size(); p++) {
        for (int m = 0; m < mol[p].size(); m++) {
            std::string type = mol[p][m].type();
            if (type == "ASN") {
                auto root = Node::create_node({p, m, "", ""});
                nodes.push_back(std::move(root));
            }
        }
    }

    for (auto& node: nodes) {
        get_value(mol, database, node);
    }

    // Remove nodes with no children.
    for (auto it = nodes.begin(); it != nodes.end();){
        if ((*it)->children.empty()) {
            it = nodes.erase(it);
        }
        else {
            ++it;
        }

    }
    for (const auto& node: nodes) {
        printTree(node.get());
    }

}
