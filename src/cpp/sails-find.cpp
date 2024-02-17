//
// Created by Jordan Dialpuri on 12/02/2024.
//

#include <optional>
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

        std::vector<double>
        operator()(const TargetFunctionBase &target_function, const std::vector<std::vector<double>> &args) const;

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
                            Sails::LinkageSet *linkage,
                            Scorer scorer,
                            float step
        ) :
                m_xmap(&xmap),
                m_linkage(linkage),
                m_scorer(scorer),
                m_step(step) {};

        int num_params() const override { return 3; }

        double operator()(const std::vector<double> &args) const override;

        LinkageSet *refine();

    private:
        clipper::Xmap<float> *m_xmap{};
        Sails::LinkageSet *m_linkage;
        Scorer m_scorer{};
        float m_step;
    };
}

float score_sugar(clipper::MMonomer &monomer, const clipper::Xmap<float> &xmap) {

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
    return score_sum / monomer.size();
}

float score_function(Sails::LinkageSet &linkage_set, clipper::Xmap<float> &map) {
    clipper::MMonomer test_monomer = linkage_set.get_donor_monomer();

    Sails::LinkageParams new_parameters = linkage_set.get_parameters();
    clipper::RTop_orth rtop = linkage_set.update_transformation(new_parameters);

    test_monomer.transform(rtop);

    float score = score_sugar(test_monomer, map);

    return score;
}

SailsFind::SailsFind(clipper::MiniMol &work_model, clipper::Xmap<float> &work_map, clipper::Xmap<float> &pred_map) {
    mol = work_model;
    pred = pred_map;
    map = work_map;
}

SailsFind::SailsFind(SailsInput &input) {
    mol = input.get_minimol();
    pred = input.get_predicted_map();
    map = input.get_work_xmap();
}

SailsFind::SailsFind(SailsInput &input, SailsOutput &output) {
    mol = input.get_minimol();
    pred = input.get_predicted_map();
    map = input.get_work_xmap();
    clipper::MiniMol mol = this->find();
    clipper::MMDBfile mf;
    mf.export_minimol(mol);
    std::string output_path = output.get_output_path();
    mf.write_file(output_path, clipper::MMDBfile::TYPE::PDB);
}

std::optional<clipper::MMonomer> SailsFind::glycosylate_protein(clipper::MMonomer &residue, clipper::MMonomer &sugar) {
    clipper::MMonomer dummy = sugar;
    Sails::ASN_NAG asnnag = {residue, sugar};

    Sails::TargetFunctionSugar refiner = {map, &asnnag, score_function, 10};
    Sails::LinkageSet *refined_linkage = refiner.refine();
    Sails::LinkageParams refined_parameters = refined_linkage->get_parameters();
    clipper::RTop_orth rtop = refined_linkage->update_transformation(refined_parameters);
    dummy.transform(rtop);

    float score = score_sugar(dummy, map);

    if (score > 0.3) {
        return dummy;
    }
    return std::nullopt;
}

std::optional<clipper::MMonomer> SailsFind::extend_sugar(clipper::MMonomer &residue, clipper::MMonomer &sugar) {
    clipper::MMonomer dummy = sugar;
    Sails::NAG_NAG nagnag = {residue, sugar};

    std::cout << "Extending sugar " << residue.type() << " => " << sugar.type() << std::endl;

    // dummy.transform(nagnag.get_transformation());
    // return dummy;

    Sails::TargetFunctionSugar refiner = {map, &nagnag, score_function, 10};

    Sails::LinkageSet *refined_linkage = refiner.refine();
    Sails::LinkageParams refined_parameters = refined_linkage->get_parameters();
    clipper::RTop_orth rtop = refined_linkage->update_transformation(refined_parameters);
    sugar.transform(rtop);

    std::cout << nagnag.get_transformation().format() << "\n"
                << nagnag.get_parameters().format() << "\n"
                << rtop.format() << "\n"
        << refined_parameters.format() << "\n"
                  <<  std::endl;

    float score = score_sugar(dummy, map);
    // std::cout << "Calculated score = " << score << std::endl;
//    if (score > 0.3) {
        return sugar;
//    }
    return std::nullopt;
}

clipper::MiniMol SailsFind::find() {

    clipper::MMonomer dummy = SailsMonomers::load_momomer("NAG");

    std::vector<GlycoSite> glycosites;

    for (int p = 0; p < mol.size(); p++) {
        for (int m = 0; m < mol[p].size(); m++) {
            std::string type = mol[p][m].type();
            if (Sails::ResidueMap.find(type) != Sails::ResidueMap.end()) {
                glycosites.emplace_back(p, m, Sails::ResidueMap[type]);
            }
        }
    }

    for (const GlycoSite &gs: glycosites) {
        switch (gs.get_residue_type()) {
            case Sails::ASN: {
                // clipper::MMonomer sugar = dummy;
                //
                // std::optional<clipper::MMonomer> positioned_sugar = glycosylate_protein(
                //         mol[gs.get_polymer_id()][gs.get_monomer_id()], sugar);
                // if (positioned_sugar.has_value()) {
                //     mol[gs.get_polymer_id()].insert(positioned_sugar.value());
                // }
                break;
            }
            case Sails::NAG: {
                clipper::MMonomer sugar = dummy;

                std::optional<clipper::MMonomer> extended_sugar = extend_sugar(
                        mol[gs.get_polymer_id()][gs.get_monomer_id()], sugar);

                if (extended_sugar.has_value()) {
                    mol[gs.get_polymer_id()].insert(extended_sugar.value());
                }
                else {
                    std::cout << "Could not build additional NAG" << std::endl;
                }
                break;
            }
            default:
                break;
        }
    }


//
//            switch (Sails::ResidueMap[type]) {
//                case Sails::ResidueType::ASN: {
//
//                    clipper::MMonomer dummy2 = dummy;
//                    Sails::ASN_NAG asnnag = {mol[p][m], dummy};
//
//                    Sails::TargetFunctionSugar refiner = {map, &asnnag, score_function,10};
//                    Sails::LinkageSet* refined_linkage = refiner.refine();
//                    Sails::LinkageParams refined_parameters = refined_linkage->get_parameters();
//                    clipper::RTop_orth rtop2 = refined_linkage->update_transformation(refined_parameters);
//                    dummy2.transform(rtop2);
//
//                    float score = score_sugar(dummy2, map);
//
//                     if (score > 0.3) {
//                         dummy2.set_id(count++);
//                         mol[p].insert(dummy2);
//                     }
//                    break;
//                }
//                case Sails::ResidueType::NAG: {
//                    clipper::MMonomer dummy2 = dummy;
//                    Sails::NAG_NAG nagnag = {mol[p][m], dummy};
//
//                }
//                default:
//                    break;
//            }
//        }
//    }
//
//    clipper::MMDBfile mf;
//    mf.export_minimol();
//    mf.write_file("debug/mol.pdb", clipper::MMDBfile::TYPE::PDB);
    return mol;
}
