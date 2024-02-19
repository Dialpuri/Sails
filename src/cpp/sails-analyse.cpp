//
// Created by Jordan Dialpuri on 19/02/2024.
//

#include "sails-analyse.h"

#include "sails-data.h"

void printTree(clipper::MiniMol& mol, const Node* root, int level = 0) {
    if(root)
    {
        std::cout << std::string(level * 4, ' ');
        std::cout << mol[root->data.get_p()].id() << "-"
        << root->data.get_mmonomer(mol).type() << "-"
        << root->data.get_mmonomer(mol).id().trim() << " / "
        << root->data.get_donor_atom() << "-" << root->data.get_acceptor_atom()
        << std::endl;

        for(const auto& child : root->children)
        {
            printTree(mol, child.get(), level + 1);
        }
    }
}

void SailsAnalyse::analyse_children(clipper::MiniMol&mol, SailsData&database, std::unique_ptr<Node>&node) {
    if (node == nullptr) {return;}

    clipper::MMonomer residue = node->data.get_mmonomer(mol);
    Sails::ResidueType type = Sails::ResidueMap[residue.type()];
    if (database.data_map.find(type) == database.data_map.end()) { std::cout << "Data on this residue is not available\n"; return;}

    Data* data = database.data_map[type];

    for (const std::string& donor: data->get_donors()) {
        int atom_index = residue.lookup(donor, clipper::MM::UNIQUE);
        if (atom_index < 0) {std::cout << "Can't find donor atom\n"; continue;}

        clipper::Coord_orth atom_coord_orth = residue[atom_index].coord_orth();
        std::vector<clipper::MAtomIndexSymmetry> atoms = ns(atom_coord_orth, 2);

        std::vector<NeighbourSearchResult> nearby_residues;
        nearby_residues.reserve(atoms.size());
        for (const auto& atom: atoms) {
            nearby_residues.emplace_back(atom.polymer(), atom.monomer(), atom.atom());
        }

        std::set<std::string> possible_links = data->get_possible_links();

        for (const auto& nearby_residue: nearby_residues) {
            if (nearby_residue.is_same(node.get())) {continue;}

            clipper::MMonomer nearby_mmonomer = nearby_residue.extract_monomer(mol);
            if (possible_links.find(nearby_mmonomer.type()) == possible_links.end()) {continue;}

            double min_dist = 1e3;
            int closest_atom_idx = -1;
            for (int a = 0; a < nearby_mmonomer.size(); a++) {
                clipper::Coord_orth co = nearby_mmonomer[a].coord_orth();
                if (double dist = clipper::Coord_orth::length(atom_coord_orth, co); dist < min_dist) {
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
        analyse_children(mol, database, child);
    }
}

SailsAnalyse::SailsAnalyse(clipper::MiniMol& mol) {

    ns = clipper::MAtomNonBond(mol, 2);
    std::vector<std::unique_ptr<Node>> nodes;

    SailsData database;

    for (int p = 0; p < mol.size(); p++) {
        for (int m = 0; m < mol[p].size(); m++) {
            if (std::string type = mol[p][m].type(); type == "ASN") {
                auto root = Node::create_node({p, m, "", ""});
                analyse_children(mol, database, root);
                if (root->children.empty()) {continue;}
                nodes.push_back(std::move(root));
            }
        }
    }

    // for (const auto& node: nodes) {
    //     printTree(mol, node.get());
    // }
    trees = nodes;
}
