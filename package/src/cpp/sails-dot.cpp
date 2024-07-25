//
// Created by Jordan Dialpuri on 24/07/2024.
//

#include "../include/sails-dot.h"



Sails::Dot::Dot(gemmi::Structure &structure) {
    JSONLoader loader = {"package/data/data.json"};
    m_database = loader.load_residue_database();
    m_structure = structure;
    m_topology = {&m_structure, m_database};
}

std::map<Sails::Glycosite, std::string> Sails::Dot::get_all_dotfiles() {
    Glycosites glycosites = find_n_glycosylation_sites(m_structure);
    std::map<Glycosite, std::string> dot_files;
    for (auto& glycosite: glycosites) {
        dot_files[glycosite] = get_dotfile(glycosite);
    }
    return dot_files;
}

std::string Sails::Dot::get_dotfile(Glycosite &glycosite) {
    Glycan glycan = m_topology.find_glycan_topology(glycosite);

    std::stringstream s;
    s << header() << std::endl;;

    for (const auto& [sugar, linked_set]: glycan.sugars) {
        s << "\"" << Utils::format_residue_from_site(sugar, &m_structure)  << "\"" <<  get_format(sugar) << std::endl;
    }
    s << footer() << std::endl;;
    for (const auto& [sugar, linked_set]: glycan.adjacency_list) {
        auto donor = Utils::format_residue_from_site(sugar->site, &m_structure);

        for (auto& linked_sugar: linked_set) {
            auto acceptor = Utils::format_residue_from_site(linked_sugar->site, &m_structure);
            s << "\"" <<  donor << "\"" << "--" << "\"" << acceptor << "\"" << std::endl;
        }
    }

    s << footer();
    return s.str();
}

std::string Sails::Dot::header() {
    return "graph Glycan {\nrankdir=RL;\npad=0.5\nnodesep=1\nrankset=2\n{";
};

std::string Sails::Dot::footer() {
    return "}";
}

// 1 [fillcolor="blue" shape=square style=filled label=""]
std::string Sails::Dot::get_format(const Glycosite &site) {
    std::string s;
    gemmi::Residue residue = Utils::get_residue_from_glycosite(site, &m_structure);

    std::string label = "";
    if (residue.name == "ASN") label = Utils::format_residue_from_site(site, &m_structure);

    s += " [";
    s += "fillcolor=\"" + m_database[residue.name].snfg_colour + "\"";
    s += "shape=\"" + m_database[residue.name].snfg_shape + "\"";
    s += "style=filled label=\"" + label + "\"]";
    return s;
}

