//
// Created by Jordan Dialpuri on 19/02/2024.
//

#ifndef SAILSANALYSE_H
#define SAILSANALYSE_H

#include <clipper/clipper-minimol.h>

#include "sails-data.h"
#include "sails-model.h"

class NodeData {
public:
    NodeData(int p, int m, const std::string& acceptor_atom, const std::string& donor_atom):
        m_p(p), m_m(m), m_acceptor_atom(acceptor_atom), m_donor_atom(donor_atom) {};

    [[nodiscard]] clipper::MMonomer get_mmonomer(clipper::MiniMol& mol) const {return mol[m_p][m_m];}
    [[nodiscard]] std::string get_acceptor_atom() const {return m_acceptor_atom;}
    [[nodiscard]] std::string get_donor_atom() const {return m_donor_atom;}
    [[nodiscard]] int get_p() const {return m_p;}
    [[nodiscard]] int get_m() const {return m_m;}
private:
    int m_p;
    int m_m;
    std::string m_acceptor_atom;
    std::string m_donor_atom;
};

struct Node {
    NodeData data;
    std::vector<std::unique_ptr<Node>> children;
    explicit Node(const NodeData& node_data): data(node_data) {}

    static std::unique_ptr<Node> create_node(const NodeData& data)
    {
        std::unique_ptr<Node> newNode = std::make_unique<Node>(data);
        newNode->data = data;
        return newNode;
    }

    static void add_child(std::unique_ptr<Node> parent, NodeData& data)
    {
        auto child = create_node(data);
        parent->children.push_back(std::move(child));
    }

    void add_child(NodeData& data)
    {
        auto child = create_node(data);
        this->children.push_back(std::move(child));
    }
};

struct NeighbourSearchResult {
    int p;
    int m;
    NeighbourSearchResult(int p_, int m_): p(p_), m(m_) {}

    bool is_same(const Node* node) const {
        return node->data.get_p() == p && node->data.get_m() == m;
    }

    clipper::MMonomer extract_monomer(clipper::MiniMol& mol) const {
        return mol[p][m];
    }

    bool operator<(const NeighbourSearchResult& rhs) const
    {
        return p < rhs.p && m < rhs.m;
    }
};

class GlycoSite {
public:
    GlycoSite(int p, int m, Sails::ResidueType& residue_type): m_polymer_id(p), m_monomer_id(m),
                                                                m_residue_type(residue_type) {}
    GlycoSite() = default;

    [[nodiscard]] Sails::ResidueType get_residue_type() const {return m_residue_type;}
    [[nodiscard]] int get_polymer_id() const {return m_polymer_id;}
    [[nodiscard]] int get_monomer_id() const {return m_monomer_id;}

private:
    int m_polymer_id;
    int m_monomer_id;
    Sails::ResidueType m_residue_type;
};


class SailsAnalyse {
public:
    bool get_value(clipper::MiniMol&mol, SailsData& database,  std::unique_ptr<Node>&node);

    explicit SailsAnalyse(clipper::MiniMol& mol);

private:
    clipper::MAtomNonBond ns;

};



#endif //SAILSANALYSE_H
