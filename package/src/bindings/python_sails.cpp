//
// Created by Jordan Dialpuri on 16/02/2024.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/map.h>


#include "../cpp/sails.cpp"
#include "../include/sails-gemmi-bindings.h"
#include "../include/sails-dot.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(sails_module, m) {
        // reflection classes
        nb::class_<Sails::NumberPair>(m, "Pair")
                        .def(nb::init<double, double>())
                        .def_rw("value1", &Sails::NumberPair::value1)
                        .def_rw("value2", &Sails::NumberPair::value2);

        nb::class_<Sails::HKL>(m, "HKL")
                        .def(nb::init<int, int, int>())
                        .def_rw("h", &Sails::HKL::h)
                        .def_rw("k", &Sails::HKL::k)
                        .def_rw("l", &Sails::HKL::l);

        nb::class_<Sails::Reflection>(m, "Reflection")
                        .def(nb::init<Sails::HKL, Sails::NumberPair>())
                        .def(nb::init<Sails::HKL, Sails::NumberPair, Sails::NumberPair>())
                        .def(nb::init<Sails::HKL, Sails::NumberPair, Sails::NumberPair, Sails::NumberPair>())
                        .def_rw("hkl", &Sails::Reflection::hkl)
                        .def_rw("f_sigf", &Sails::Reflection::f_sigf)
                        .def_rw("fwt_phwt", &Sails::Reflection::fwt_phwt)
                        .def_rw("delfwt_phdelwt", &Sails::Reflection::delfwt_phdelwt);


        // gemmi Structure
        nb::class_<gemmi::Structure>(m, "Structure")
                        .def(nb::init<>())
                        .def_rw("models", &gemmi::Structure::models)
                        .def("cell", [](const gemmi::Structure &structure) {
                                return Sails::Cell(structure.cell);
                        })
                        .def("set_cell", [](gemmi::Structure &structure, const Sails::Cell &cell) {
                                structure.cell = gemmi::UnitCell(cell.a, cell.b, cell.c, cell.alpha, cell.beta,
                                                                 cell.gamma);
                        });

        nb::bind_vector<std::vector<gemmi::Model> >(m, "Models");

        nb::class_<gemmi::Model>(m, "Model")
                        .def(nb::init<>())
                        .def_rw("name", &gemmi::Model::name)
                        .def_rw("chains", &gemmi::Model::chains);
        nb::bind_vector<std::vector<gemmi::Chain> >(m, "Chains");

        nb::class_<gemmi::Chain>(m, "Chain")
                        .def(nb::init<>())
                        .def_rw("residues", &gemmi::Chain::residues)
                        .def_rw("name", &gemmi::Chain::name);
        nb::bind_vector<std::vector<gemmi::Residue> >(m, "Residues");

        nb::class_<gemmi::Residue>(m, "Residue")
                        .def(nb::init<>())
                        .def_rw("atoms", &gemmi::Residue::atoms)
                        .def_rw("name", &gemmi::Residue::name)
                        .def_rw("seqid", &gemmi::Residue::seqid)
                        .def_rw("subchain", &gemmi::Residue::subchain)
                        .def("get_label_seq", [ ](gemmi::Residue &residue) {
                                return residue.label_seq.value;
                        })
                        .def("set_label_seq", [ ](gemmi::Residue &residue, const int label_seq) {
                                residue.label_seq = gemmi::SeqId::OptionalNum(label_seq);
                        })
                        .def_rw("entity_id", &gemmi::Residue::entity_id);
        nb::bind_vector<std::vector<gemmi::Atom> >(m, "Atoms");

        nb::class_<gemmi::Atom>(m, "Atom")
                        .def(nb::init<>())
                        .def_rw("pos", &gemmi::Atom::pos)
                        .def_rw("b_iso", &gemmi::Atom::b_iso)
                        .def_rw("altloc", &gemmi::Atom::altloc)
                        .def("set_altloc", [](gemmi::Atom& atom, std::string& altloc) {
                            atom.altloc = altloc[0];
                        })
                        .def_rw("element", &gemmi::Atom::element)
                        .def_rw("occ", &gemmi::Atom::occ)
                        .def_rw("name", &gemmi::Atom::name);

        nb::class_<gemmi::Position>(m, "Position")
                        .def(nb::init<double, double, double>())
                        .def_ro("x", &gemmi::Position::x)
                        .def_ro("y", &gemmi::Position::y)
                        .def_ro("z", &gemmi::Position::z);

        nb::class_<gemmi::Element>(m, "Element")
                        .def(nb::init<std::string &>())
                        .def_prop_ro("name", &gemmi::Element::name);

        nb::class_<gemmi::SeqId>(m, "SeqId")
                        .def(nb::init<int, char>())
                        .def("num", [](gemmi::SeqId &s) {
                                return s.num.value;
                        })
                        .def("icode", [](gemmi::SeqId &s) {
                                return s.icode;
                        });

        nb::bind_vector<std::vector<Sails::Reflection> >(m, "Reflections");

        nb::class_<Sails::Cell>(m, "Cell")
                        .def(nb::init<double, double, double, double, double, double>())
                        .def_ro("a", &Sails::Cell::a)
                        .def_ro("b", &Sails::Cell::b)
                        .def_ro("c", &Sails::Cell::c)
                        .def_ro("alpha", &Sails::Cell::alpha)
                        .def_ro("beta", &Sails::Cell::beta)
                        .def_ro("gamma", &Sails::Cell::gamma);

        nb::class_<Sails::MTZ>(m, "MTZ")
                        .def(nb::init<std::vector<Sails::Reflection> &, Sails::Cell &, std::string &>())
                        .def_ro("reflections", &Sails::MTZ::reflections)
                        .def_ro("cell", &Sails::MTZ::cell)
                        .def_ro("spacegroup", &Sails::MTZ::spacegroup);

        nb::class_<Sails::Output>(m, "SailsOutput")
                        .def_ro("structure", &Sails::Output::structure)
                        .def_ro("mtz", &Sails::Output::mtz)
                        .def_ro("log", &Sails::Output::log)
                        .def_ro("snfgs", &Sails::Output::snfgs);

        nb::class_<Sails::Glycosite>(m, "GlycoSite")
                        .def(nb::init<int, int, int>())
                        .def_rw("model_idx", &Sails::Glycosite::model_idx)
                        .def_rw("chain_idx", &Sails::Glycosite::chain_idx)
                        .def_rw("residue_idx", &Sails::Glycosite::residue_idx)
                        .def_rw("atom_idx", &Sails::Glycosite::atom_idx);

        nb::class_<Sails::Dot>(m, "Dot")
                        .def(nb::init<gemmi::Structure &>())
                        .def("get_all_dotfiles", &Sails::Dot::get_all_dotfiles)
                        .def("get_dotfile", &Sails::Dot::get_dotfile);

        nb::enum_<gemmi::AxisOrder>(m, "AxisOrder")
                .value("XYZ", gemmi::AxisOrder::XYZ)
                .value("ZYX", gemmi::AxisOrder::ZYX)
                .value("Unknown", gemmi::AxisOrder::Unknown);

    nb::class_<gemmi::Grid<float>>(m, "Grid")
                        .def(nb::init<>())
                        .def("set_data",  [](gemmi::Grid<float>& grid, nb::ndarray<>& arr) {
                            grid.data = {static_cast<float *>(arr.data()), static_cast<float *>(arr.data())+arr.size()};
                        })
                        .def_rw("nu", &gemmi::Grid<>::nu)
                        .def_rw("nv", &gemmi::Grid<>::nv)
                        .def_rw("nw", &gemmi::Grid<>::nw)
                        .def_rw("cell", &gemmi::Grid<>::unit_cell)
                        .def("set_spacegroup", [](gemmi::Grid<float>& grid, std::string& spacegroup) {
                            grid.spacegroup = gemmi::find_spacegroup_by_name(spacegroup);
                        })
                        .def("set_cell", [](gemmi::Grid<float>& grid, Sails::Cell& cell) {
                            grid.set_unit_cell(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
                        })
                        .def("set_axis_order", [](gemmi::Grid<float>& grid, gemmi::AxisOrder axis_order) {
                            grid.axis_order = axis_order;
                        })
                        .def("set_spacing", [](gemmi::Grid<float>& grid, double xspacing, double yspacing, double zspacing) {
                            grid.spacing[0] = xspacing;
                            grid.spacing[1] = yspacing;
                            grid.spacing[2] = zspacing;
                        })

                        ;


    m.def("n_glycosylate",
          nb::overload_cast<gemmi::Structure &, Sails::MTZ &, int, std::string &, bool>(&n_glycosylate), "structure"_a,
          "mtz"_a, "cycles"_a, "resource_dir"_a, "verbose"_a);
    m.def("c_glycosylate",
          nb::overload_cast<gemmi::Structure &, Sails::MTZ &, int, std::string &, bool>(&c_glycosylate), "structure"_a,
          "mtz"_a, "cycles"_a, "resource_dir"_a, "verbose"_a);
    m.def("o_mannosylate",
          nb::overload_cast<gemmi::Structure &, Sails::MTZ &, int, std::string &, bool>(&o_mannosylate), "structure"_a,
          "mtz"_a, "cycles"_a, "resource_dir"_a, "verbose"_a);

    m.def("n_glycosylate",
          nb::overload_cast<gemmi::Structure &, gemmi::Grid<> &, int, std::string &, bool>(&n_glycosylate),
          "structure"_a, "grid"_a, "cycles"_a, "resource_dir"_a, "verbose"_a);
    m.def("c_glycosylate",
          nb::overload_cast<gemmi::Structure &, gemmi::Grid<> &, int, std::string &, bool>(&c_glycosylate),
          "structure"_a, "grid"_a, "cycles"_a, "resource_dir"_a, "verbose"_a);
    m.def("o_mannosylate",
          nb::overload_cast<gemmi::Structure &, gemmi::Grid<> &, int, std::string &, bool>(&o_mannosylate),
          "structure"_a, "grid"_a, "cycles"_a, "resource_dir"_a, "verbose"_a);

    m.def("find_all_wurcs", &find_all_wurcs, "structure"_a, "resource_dir"_a);
    m.def("find_wurcs", &find_wurcs, "structure"_a, "chain"_a, "seqid"_a,  "resource_dir"_a);
    m.def("model_wurcs", &model_wurcs, "structure"_a, "wurcs"_a, "chain"_a, "seqid"_a, "resource_dir"_a);

    m.def("morph", &morph, "structure"_a, "wurcs"_a, "chain"_a, "seqid"_a, "resource_dir"_a);

    m.def("test_snfg", &test);

    m.def("get_snfg", &get_snfg, "chain"_a, "seqid"_a , "structure"_a, "resource_dir"_a);
    m.def("get_all_snfgs", &get_all_snfgs, "structure"_a, "resource_dir"_a);
}
