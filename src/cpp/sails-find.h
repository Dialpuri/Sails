//
// Created by Jordan Dialpuri on 12/02/2024.
//
#pragma once

#ifndef SAILS_FIND_H
#define SAILS_FIND_H

#include <clipper/clipper-minimol.h>
#include <optional>
#include "sails-lib.h"
#include "sails-model.h"
#include "sails-util.h"

class GlycoSite {
public:
    GlycoSite(int p, int m, Sails::ResidueType& residue_type): m_polymer_id(p), m_monomer_id(m),
                                                                m_residue_type(residue_type) {}
    GlycoSite() = default;

    Sails::ResidueType get_residue_type() const {return m_residue_type;}
    [[nodiscard]] int get_polymer_id() const {return m_polymer_id;}
    [[nodiscard]] int get_monomer_id() const {return m_monomer_id;}

private:
    int m_polymer_id;
    int m_monomer_id;
    Sails::ResidueType m_residue_type;
};

/**
*
*/
class SailsFind {
public:
    SailsFind(clipper::MiniMol& work_model, clipper::Xmap<float>& work_map, clipper::Xmap<float>& pred_map);
    SailsFind(SailsInput& input);
    SailsFind(SailsInput& input, SailsOutput& output);

    clipper::MiniMol find();

private:
    std::optional<clipper::MMonomer> glycosylate_protein(clipper::MMonomer &residue, clipper::MMonomer &sugar);
    std::optional<clipper::MMonomer> extend_sugar(clipper::MMonomer &residue, clipper::MMonomer &sugar);

private:
    clipper::MiniMol mol;
    clipper::Xmap<float> pred;
    clipper::Xmap<float> map;

};
#endif //SAILS_FIND_H
