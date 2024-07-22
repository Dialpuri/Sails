//
// Created by Jordan Dialpuri on 22/07/2024.
//

#include "../include/sails-cif.h"

std::vector<Sails::LinkRecord>  Sails::generate_link_records(gemmi::Structure *structure, Glycosites *glycosites,
                                                             Topology *topology) {
    std::vector<LinkRecord> links;
    int link_id = 0;
    for (auto &glycosite: *glycosites) {
        Glycan glycan = topology->find_glycan_topology(glycosite);
        std::vector<Linkage> list = glycan.linkage_list;

        for (const auto &l: list) {
            std::string id = "covalent" + std::to_string(link_id++);
            gemmi::Chain c1 = Utils::get_chain_from_glycosite(l.donor_sugar->site, structure);
            gemmi::Chain c2 = Utils::get_chain_from_glycosite(l.acceptor_sugar->site, structure);

            gemmi::Residue r1 = Utils::get_residue_from_glycosite(l.donor_sugar->site, structure);
            gemmi::Residue r2 = Utils::get_residue_from_glycosite(l.acceptor_sugar->site, structure);

            gemmi::Atom *a1 = &r1.get(l.donor_atom)[0];
            gemmi::Atom *a2 = &r2.get(l.acceptor_atom)[0];

            LinkRecord link = {id, c1, c2, r1, r2, *a1, *a2};
            links.push_back(link);
        }
    }
    return links;
}
