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

void Sails::add_links_to_structure(gemmi::Structure *structure, std::vector<Sails::LinkRecord> &link_records) {

    std::vector<gemmi::Connection> connections;

    for (auto& link_record : link_records) {
        gemmi::Connection connection;
        connection.type = gemmi::Connection::Covale;
        gemmi::AtomAddress a1;
        a1.chain_name = link_record.chain1.name;
        gemmi::ResidueId resid1;
        resid1.name = link_record.residue1.name;
        resid1.seqid = link_record.residue1.seqid;
        a1.res_id = resid1;
        a1.atom_name = link_record.atom1.name;

        gemmi::AtomAddress a2;
        a2.chain_name = link_record.chain2.name;
        gemmi::ResidueId resid2;
        resid2.name = link_record.residue2.name;
        resid2.seqid = link_record.residue2.seqid;
        a2.res_id = resid2;
        a2.atom_name = link_record.atom2.name;

        double distance = (link_record.atom1.pos - link_record.atom2.pos).length();

        connection.partner1 = a1;
        connection.partner2 = a2;
        connection.reported_distance = distance;
        connections.push_back(connection);

   }

    structure->connections = connections;
}
