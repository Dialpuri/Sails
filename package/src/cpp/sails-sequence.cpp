//
// Created by Jordan Dialpuri on 06/07/2024.
//
#include "src/include/sails-sequence.h"

Sails::Glycosites Sails::find_n_glycosylation_sites(const gemmi::Structure &structure) {

    std::vector<Glycosite> glycosites;
    auto models = structure.children();
    for (int model_idx = 0; model_idx < models.size(); model_idx++) {

        auto chains = models[model_idx].children();
        for (int chain_idx = 0; chain_idx < chains.size(); chain_idx++) {

            auto residues = chains[chain_idx].children();
            for (int residue_idx = 0; residue_idx < residues.size() - 2; residue_idx++) {

                char first = gemmi::find_tabulated_residue(residues[residue_idx].name).one_letter_code;
                if (first != 'N') {continue;}

                char third = gemmi::find_tabulated_residue(residues[residue_idx + 2].name).one_letter_code;
                if (third != 'S' && third != 'T') {continue;}

                char second = gemmi::find_tabulated_residue(residues[residue_idx + 1].name).one_letter_code;
                if (second == 'P') {continue;}

                glycosites.emplace_back(model_idx, chain_idx, residue_idx);
            }
        }
    }

    return glycosites;
}

Sails::Glycosites Sails::find_c_glycosylation_sites(const gemmi::Structure &structure) {
    std::vector<Glycosite> glycosites;
    auto models = structure.children();
    for (int model_idx = 0; model_idx < models.size(); model_idx++) {

        auto chains = models[model_idx].children();
        for (int chain_idx = 0; chain_idx < chains.size(); chain_idx++) {

            auto residues = chains[chain_idx].children();
            for (int residue_idx = 0; residue_idx < residues.size() - 3; residue_idx++) {

                char first = gemmi::find_tabulated_residue(residues[residue_idx].name).one_letter_code;
                if (first != 'W') {continue;}

                char fourth = gemmi::find_tabulated_residue(residues[residue_idx + 3].name).one_letter_code;
                if (fourth != 'W') {continue;}

                glycosites.emplace_back(model_idx, chain_idx, residue_idx);
                glycosites.emplace_back(model_idx, chain_idx, residue_idx+3);
            }
        }
    }

    return glycosites;
}
