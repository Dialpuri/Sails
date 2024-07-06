//
// Created by Jordan Dialpuri on 06/07/2024.
//

#include "../include/sails-topology.h"

std::optional<Sails::Glycan>
Sails::find_glycan_topology(gemmi::Structure &structure, const Sails::Glycosite &glycosite) {

    constexpr double search_radius = 2;
    gemmi::NeighborSearch neighbor_search = {structure.models[glycosite.model_idx], structure.cell, search_radius};
    neighbor_search.populate(false);

    gemmi::Residue amino_acid = structure.models[glycosite.model_idx].chains[glycosite.chain_idx].residues[glycosite.residue_idx];

    auto near_atoms = neighbor_search.find_atoms(amino_acid.sole_atom("ND2").pos, '\0', 0.0, search_radius);
    std::cout << amino_acid.seqid.num.value << " has " << near_atoms.size() << " nearby" << std::endl;

    return std::nullopt;
}
