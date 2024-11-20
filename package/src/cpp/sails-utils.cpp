//
// Created by Jordan Dialpuri on 08/07/2024.
//

#include <utility>

#include "../include/sails-utils.h"


std::string Sails::Utils::get_environment_variable(std::string const &key) {
    char *val = getenv(key.c_str());
    return val == nullptr ? std::string("") : std::string(val);
}

void Sails::Utils::print_vector(const gemmi::Vec3 &vec) {
    std::cout << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
}

void Sails::Utils::print_vector(const gemmi::Vec3 &vec, const std::string &name) {
    std::cout << name << ": (" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
}

void Sails::Utils::print_position(const gemmi::Position &position) {
    std::cout << "(" << position.x << ", " << position.y << ", " << position.z << ")" << std::endl;
}

std::string Sails::Utils::format_residue_key(const gemmi::Residue *residue) {
    return residue->name + "-" + residue->seqid.str();
}

std::string Sails::Utils::format_site_key(const Glycosite &site) {
    return std::to_string(site.model_idx) + "-" +
           std::to_string(site.chain_idx) + "-" +
           std::to_string(site.residue_idx) + "-" +
           std::to_string(site.atom_idx);
}

std::string Sails::Utils::format_residue_from_site(const Glycosite &glycosite, gemmi::Structure *structure) {
    gemmi::Residue* residue = get_residue_ptr_from_glycosite(glycosite, structure);
    gemmi::Chain* chain = get_chain_ptr_from_glycosite(glycosite, structure);
    return chain->name + "-" + residue->name + "-" + residue->seqid.num.str();

}

gemmi::Chain Sails::Utils::get_chain_from_glycosite(const Glycosite &site, const gemmi::Structure* structure) {
    return structure->models[site.model_idx].chains[site.chain_idx];
}

gemmi::Residue Sails::Utils::get_residue_from_glycosite(const Glycosite &site, const gemmi::Structure *structure) {
    return structure->models[site.model_idx].chains[site.chain_idx].residues[site.residue_idx];
}

gemmi::Residue *Sails::Utils::
get_residue_ptr_from_glycosite(const Glycosite &site, gemmi::Structure *structure) {
    return &structure->models[site.model_idx].chains[site.chain_idx].residues[site.residue_idx];
}

gemmi::Chain *Sails::Utils::
get_chain_ptr_from_glycosite(const Glycosite &site, gemmi::Structure *structure) {
    return &structure->models[site.model_idx].chains[site.chain_idx];
}

gemmi::Atom Sails::Utils::get_atom_from_glycosite(const Glycosite &site, const gemmi::Structure* structure) {
    if (site.atom_idx == -1) { throw std::runtime_error("Site has not been initialised from a Mark"); }
    return structure->models[site.model_idx].chains[site.chain_idx].residues[site.residue_idx].atoms[site.atom_idx];
}

std::string Sails::Utils::linkage_to_id(const Sails::LinkageData &data) {
    return data.donor + "-" + std::to_string(data.donor_number) + "," + std::to_string(data.acceptor_number) + "-" +
           data.acceptor;
}

void Sails::Utils::save_residues_to_file(std::vector<gemmi::Residue> residues, const std::string &path) {
    gemmi::Structure structure;
    gemmi::Model model = gemmi::Model("A");
    gemmi::Chain chain = gemmi::Chain("A");
    chain.append_residues(std::move(residues));
    model.chains = {chain};
    structure.models = {model};

    std::ofstream file;
    file.open(path);
    gemmi::write_pdb(structure, file);
    file.close();
}

void Sails::Utils::save_grid_to_file(const gemmi::Grid<> &grid, const std::string &path) {
    gemmi::Ccp4<> map;
    map.grid = grid;
    map.update_ccp4_header();
    map.write_ccp4_map(path);
}

void Sails::Utils::save_structure_to_file(const gemmi::Structure &structure, const std::string &path) {
    std::ofstream of;
    of.open(path);
    write_cif_to_stream(of, make_mmcif_document(structure));
    of.close();
}

std::string Sails::Utils::get_next_string(const std::string& input) {
    std::string result = input;
    int i = result.length() - 1;

    while (i >= 0) {
        if (result[i] == 'Z') {
            result[i] = 'A';
            i--;
        } else {
            result[i] += 1;
            return result;
        }
    }
    return "A" + result;
}

bool Sails::Utils::file_exists(const std::string &path) {
    std::ifstream f;
    f.open(path);
    const bool good = f.good();
    f.close();
    return good;
}

std::vector<std::string> Sails::Utils::split(const std::string &string, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream token_stream(string);
    while (std::getline(token_stream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
