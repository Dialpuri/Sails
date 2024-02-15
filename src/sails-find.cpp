//
// Created by Jordan Dialpuri on 12/02/2024.
//

#include "sails-find.h"

namespace Sails {

    class TargetFunctionBase {
    public:
        TargetFunctionBase() {}

        virtual ~TargetFunctionBase() {}

        virtual int num_params() const = 0;

        virtual double operator()(const std::vector<double> &args) const = 0;
    };

    class SimplexOptimiser {
    public:
        enum TYPE {
            NORMAL, GRADIENT
        };

        SimplexOptimiser(double tolerance = 1e-3, int max_cycles = 50, TYPE type = NORMAL, bool debug = false) :
                m_tolerance(tolerance), m_max_cycles(max_cycles), m_type(type), m_debug_mode(debug) {};

        std::vector<double> operator()(const TargetFunctionBase &target_function, const std::vector<std::vector<double>> &args) const;

        void debug() { m_debug_mode = true; }

    private:
        double m_tolerance;
        double m_max_cycles;
        TYPE m_type;
        bool m_debug_mode;
    };

    class TargetFunctionSugar : public TargetFunctionBase {
    public:
        typedef float (*Scorer)(Sails::LinkageSet &, clipper::Xmap<float> &);

        TargetFunctionSugar() = default;

        TargetFunctionSugar(clipper::Xmap<float> &xmap,
                            Sails::LinkageSet* linkage,
                            Scorer scorer,
                            float step
        ) :
                m_xmap(&xmap),
                m_linkage(linkage),
                m_scorer(scorer),
                m_step(step) {};

        int num_params() const override { return 3; }

        double operator()(const std::vector<double> &args) const override;

        LinkageSet * refine();

    private:
        clipper::Xmap<float> *m_xmap{};
        Sails::LinkageSet* m_linkage;
        Scorer m_scorer{};
        float m_step;
    };
}

float score_sugar(clipper::MMonomer& monomer, const clipper::Xmap<float>&xmap) {

    if (monomer.size() == 0) {
        std::cout << "No atoms in the positioned monomer" << std::endl;
        return -10;
    }
    float score_sum = 0;
    for (int a = 0; a < monomer.size(); a++) {
        clipper::Coord_orth coordinate = monomer[a].coord_orth();
        clipper::Coord_frac coordinate_frac = coordinate.coord_frac(xmap.cell());
        score_sum += xmap.interp<clipper::Interp_cubic>(coordinate_frac);
    }
    return score_sum/monomer.size();
}

float score_function(Sails::LinkageSet& linkage_set, clipper::Xmap<float>& map) {
    clipper::MMonomer test_monomer = linkage_set.get_donor_monomer();

    Sails::LinkageParams new_parameters = linkage_set.get_parameters();
    clipper::RTop_orth rtop = linkage_set.update_transformation(new_parameters);

    test_monomer.transform(rtop);

    float score = score_sugar(test_monomer, map);

    return score;
}

SailsFind::SailsFind(clipper::MiniMol& work_model, clipper::Xmap<float>& work_map, clipper::Xmap<float>& pred_map) {
    mol = work_model;
    pred = pred_map;
    map = work_map;
}

clipper::MiniMol SailsFind::find() {

    clipper::MMonomer dummy = SailsMonomers::load_momomer("NAG");
    std::cout << dummy.size() << std::endl;
//    clipper::MiniMol dbg = {mol.spacegroup(), mol.cell()};
//    clipper::MPolymer mp;
//    mp.set_id("A");
//    clipper::MMonomer mm;
//    mm.set_id("1");
//    mm.set_type("DUM");

    int count = 0;
    for (int p = 0; p < mol.size(); p++) {
        for (int m = 0; m < mol[p].size(); m++) {
            std::string type = mol[p][m].type();
            if (Sails::AminoAcidMap.find(type) == Sails::AminoAcidMap.end()) {
                continue;
            }

            switch (Sails::AminoAcidMap[type]) {
                case Sails::AminoAcidType::ASN: {

                    clipper::MMonomer dummy2 = dummy;
                    Sails::ASN_NAG asnnag = {mol[p][m], dummy};
//                    clipper::RTop_orth rtop = asnnag.get_transformation();
//                    dummy.transform(rtop);
//                    mp.insert(dummy);
//                    std::cout << "Original Parameters = " << asnnag.get_parameters().format() << std::endl;

                    Sails::TargetFunctionSugar refiner = {map, &asnnag, score_function,10};
                    Sails::LinkageSet* refined_linkage = refiner.refine();
                    Sails::LinkageParams refined_parameters = refined_linkage->get_parameters();
                    clipper::RTop_orth rtop2 = refined_linkage->update_transformation(refined_parameters);
                    dummy2.transform(rtop2);
//                    mp.insert(dummy2);

//                    std::cout << "Final Parameters = " << refined_parameters.format() << std::endl;

                    float score = score_sugar(dummy2, map);

                     if (score > 0.3) {
                         dummy2.set_id(count++);
                         mol[p].insert(dummy2);
                     }

//                    mp.insert( mol[p][m]);
//                     dbg.insert(mp);
//                     clipper::MMDBfile mf;
//                     mf.export_minimol(dbg);
//                     mf.write_file("debug/mol.pdb");
//                    exit(-1);
                    break;
                }
                default:
                    break;
            }
        }
    }

    clipper::MMDBfile mf;
    mf.export_minimol(mol);
    mf.write_file("debug/mol.pdb", clipper::MMDBfile::TYPE::PDB);
    return clipper::MiniMol();
}
