//
// Created by Jordan Dialpuri on 12/02/2024.
//

#include "sails-find.h"


SailsFind::SailsFind(clipper::MiniMol& work_model) {
    mol = work_model;
}

clipper::MiniMol SailsFind::find() {

    clipper::MMonomer dummy = SailsMonomers::load_momomer("NAG");

    clipper::MiniMol dbg;
    clipper::MPolymer mp;
    mp.set_id("A");
    clipper::MMonomer mm;
    mm.set_id("1");
    mm.set_type("DUM");

    for (int p = 0; p < mol.size(); p++) {
        for (int m = 0; m < mol[p].size(); m++) {
            std::string type = mol[p][m].type();
            if (type == "ASN") {
                ASN asn = ASN(mol[p][m]);

                clipper::Coord_orth c1 = asn.get_position1();
                clipper::Coord_orth o5 = asn.get_position2();

                NAG nag = NAG(dummy, c1, o5);
                dummy = nag.get_positioned_monomer();

                // mm.insert(SailsUtil::create_atom(c1, "1"));
                // mm.insert(SailsUtil::create_atom(o5, "2"));

                mp.insert(dummy);
                dbg.insert(mp);
                clipper::MMDBfile mf;
                mf.export_minimol(dbg);
                mf.write_file("debug/mol.pdb");

                exit(1);
            }

        }
    }

    return clipper::MiniMol();
}
