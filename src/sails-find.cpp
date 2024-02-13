//
// Created by Jordan Dialpuri on 12/02/2024.
//

#include "sails-find.h"


SailsFind::SailsFind(clipper::MiniMol& work_model, clipper::Xmap<float>& work_map, clipper::Xmap<float>& pred_map) {
    mol = work_model;
    pred = pred_map;
    map = work_map;
}

clipper::MiniMol SailsFind::find() {

    clipper::MMonomer dummy = SailsMonomers::load_momomer("DUM");

    clipper::MiniMol dbg;
    clipper::MPolymer mp;
    mp.set_id("A");
    clipper::MMonomer mm;
    mm.set_id("1");
    mm.set_type("DUM");

    int count = 0;
    for (int p = 0; p < mol.size(); p++) {
        for (int m = 0; m < mol[p].size(); m++) {
            std::string type = mol[p][m].type();
            if (Sails::AminoAcidMap.find(type) == Sails::AminoAcidMap.end()) {
                continue;
            }

            switch (Sails::AminoAcidMap[type]) {
                case Sails::AminoAcidType::ASN: {

                    Sails::ASN_NAG asnnag = {mol[p][m], dummy};
                    clipper::RTop_orth rtop = asnnag.get_transformation();
                    dummy.transform(rtop);

                    // Sails::Asparagine asn = Sails::Asparagine(mol[p][m]);
                    //
                    // clipper::Coord_orth c1 = asn.get_position1();
                    // clipper::Coord_orth o5 = asn.get_position2();
                    // clipper::Coord_orth c5 = asn.get_position3();
                    //
                    // Sails::NAG nag = Sails::NAG(dummy, c1, o5, c5);
                    // float score = nag.score_sugar(map);
                    // if (score > 0.3) {
                    //     dummy = nag.get_positioned_monomer();
                    //     dummy.set_id(count++);
                    mol[p].insert(dummy);
                    // }

                    // mp.insert(dummy);
                    // dbg.insert(mp);
                    // clipper::MMDBfile mf;
                    // mf.export_minimol(dbg);
                    // mf.write_file("debug/mol.pdb");

                    break;
                }
                default:
                    break;
            }
        }
    }

    clipper::MMDBfile mf;
    mf.export_minimol(mol);
    mf.write_file("debug/mol.cif", clipper::MMDBfile::TYPE::CIF);
    return clipper::MiniMol();
}
