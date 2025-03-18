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
            std::string id = "covale" + std::to_string(link_id++);
            gemmi::Chain c1 = Utils::get_chain_from_glycosite(l.donor_sugar->site, structure);
            gemmi::Chain c2 = Utils::get_chain_from_glycosite(l.acceptor_sugar->site, structure);

            gemmi::Residue r1 = Utils::get_residue_from_glycosite(l.donor_sugar->site, structure);
            gemmi::Residue r2 = Utils::get_residue_from_glycosite(l.acceptor_sugar->site, structure);

            gemmi::Atom *a1 = &r1.get(l.donor_atom)[0];
            gemmi::Atom *a2 = &r2.get(l.acceptor_atom)[0];

            const float distance = (a1->pos - a2->pos).length();

            LinkRecord link = {id, c1, c2, r1, r2, *a1, *a2, distance};
            links.push_back(link);
        }
    }
    return links;
}


void Sails::add_link_records_to_structure(gemmi::Structure *structure, std::vector<Sails::LinkRecord> &link_records) {
    for (auto &link: link_records) {
        gemmi::AtomAddress a1 = {link.chain1.name, link.residue1.seqid, link.residue1.name, link.atom1.name};
        gemmi::AtomAddress a2 = {link.chain2.name, link.residue2.seqid, link.residue2.name, link.atom2.name};
        gemmi::Connection connection = {
            link.id,
            "",
            gemmi::Connection::Covale,
            gemmi::Asu::Any,
            a1,
            a2,
            link.distance
        };
        structure->connections.push_back(connection);
    }
}
