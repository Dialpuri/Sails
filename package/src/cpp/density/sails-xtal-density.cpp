//
// Created by Jordan Dialpuri on 13/08/2024.
//

#include "../../include/density/sails-density.h"
#include "../../include/density/sails-xtal-density.h"

Sails::XtalDensity::XtalDensity(gemmi::Mtz &mtz) {
    m_mtz = std::move(mtz);
    initialise_hkl();
    load_hkl("FP", "SIGFP");
}

Sails::XtalDensity::XtalDensity(gemmi::Mtz &mtz, const std::string& F, const std::string& SIGF) {
    m_mtz = std::move(mtz);
    initialise_hkl();
    load_hkl(F, SIGF);
}


void Sails::XtalDensity::initialise_hkl() {
    m_resolution = clipper::Resolution(m_mtz.resolution_high());
    const gemmi::UnitCell *cell = &m_mtz.cell;
    m_cell = clipper::Cell(clipper::Cell_descr(
            cell->a, cell->b, cell->c, cell->alpha, cell->beta, cell->gamma
    ));
    m_spacegroup = clipper::Spacegroup(clipper::Spgr_descr(m_mtz.spacegroup_name));
    m_hkl_info = {m_spacegroup, m_cell, m_resolution, true};
    m_grid_sampling = {m_spacegroup, m_cell, m_resolution};
}

void Sails::XtalDensity::load_hkl(const std::string &f, const std::string &sig_f) {
    gemmi::AsuData<gemmi::ValueSigma<float> > values = gemmi::make_asu_data<gemmi::ValueSigma<float>, 2>(
            m_mtz, {f, sig_f}, false);

    std::vector<clipper::HKL> hkls;
    hkls.reserve(values.v.size());
    std::map<std::vector<int>, std::pair<float, float> > hkl_map = {};
    for (auto &[hkl, value]: values.v) {
        hkls.emplace_back(hkl[0], hkl[1], hkl[2]);
        hkl_map.insert({
                               {hkl[0], hkl[1], hkl[2]},
                               std::make_pair(value.value, value.sigma)
                       });
    }

    m_hkl_info.add_hkl_list(hkls);

    m_fobs = clipper::HKL_data<clipper::data32::F_sigF>(m_hkl_info);

    for (HRI ih = m_fobs.first(); !ih.last(); ih.next()) {
        clipper::HKL hkl = ih.hkl();
        const auto key = {hkl.h(), hkl.k(), hkl.l()};
        auto [f, sigf] = hkl_map[key];
        const clipper::data32::F_sigF fphi = {f, sigf};
        m_fobs[ih] = fphi;
    }
}

gemmi::Grid<> Sails::XtalDensity::load_grid(const gemmi::Mtz &mtz, const std::string &f_col, const std::string &phi_col,
                                        bool normalise) {
    constexpr std::array<int, 3> null_size = {0, 0, 0};
    constexpr double sample_rate = 0;
    constexpr auto order = gemmi::AxisOrder::XYZ;

    const gemmi::Mtz::Column &f = mtz.get_column_with_label(f_col);
    const gemmi::Mtz::Column &phi = mtz.get_column_with_label(phi_col);
    const gemmi::FPhiProxy fphi(gemmi::MtzDataProxy{mtz}, f.idx, phi.idx);
    gemmi::Grid<> grid = gemmi::transform_f_phi_to_map2<float>(fphi, null_size, sample_rate, null_size, order);
    if (normalise) grid.normalize();

    return grid;
}


void Sails::XtalDensity::form_atom_list(const gemmi::Structure &structure, std::vector<clipper::Atom> &atoms) {
    for (const auto &model: structure.models) {
        for (const auto &chain: model.chains) {
            for (const auto &residue: chain.residues) {
                for (const auto &atom: residue.atoms) {
                    gemmi::Position pos = atom.pos;
                    clipper::Atom clipper_atom;
                    clipper_atom.set_element(clipper::String(atom.element.name()));
                    clipper_atom.set_coord_orth(clipper::Coord_orth(pos.x, pos.y, pos.z));
                    clipper_atom.set_occupancy(atom.occ);
                    clipper_atom.set_u_iso(clipper::Util::b2u(atom.b_iso));
                    atoms.emplace_back(clipper::MAtom(clipper_atom));
                }
            }
        }
    }
}


void Sails::XtalDensity::recalculate_map(gemmi::Structure &structure) {
    std::vector<clipper::Atom> atoms;
    form_atom_list(structure, atoms);

    clipper::Xmap<float> calculated_map = {m_spacegroup, m_cell, m_grid_sampling};
    clipper::EDcalc_iso<float> ed_calc = {m_resolution.limit()};
    clipper::HKL_data<clipper::data32::F_phi> fphic(m_hkl_info);

    ed_calc(calculated_map, atoms);
    calculated_map.fft_to(fphic);

    clipper::HKL_data<clipper::data32::Flag> modeflag(m_fobs);
    for (HRI ih = modeflag.first(); !ih.last(); ih.next())
        if (!m_fobs[ih].missing())
            modeflag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
        else
            modeflag[ih].flag() = clipper::SFweight_spline<float>::NONE;

    clipper::HKL_data<clipper::data32::F_phi> fdiff(m_fobs);
    clipper::HKL_data<clipper::data32::F_phi> fbest(m_fobs);
    clipper::HKL_data<clipper::data32::Phi_fom> phiw(m_fobs);

    clipper::SFweight_spline<float> sfw(m_hkl_info.num_reflections(), 20);
    bool success = sfw(fbest, fdiff, phiw, m_fobs, fphic, modeflag);
    if (!success) throw std::runtime_error("Sigma-A calculation failed");

    m_best_map = {m_spacegroup, m_cell, m_grid_sampling};
    m_best_map.fft_from(fbest);

    std::vector<float> recalculated_data;
    int num_columns = 6;
    recalculated_data.reserve(fdiff.data_size() * num_columns);
    for (HRI ih = m_fobs.first(); !ih.last(); ih.next()) {
        clipper::HKL hkl = ih.hkl();
        if (hkl.h() == 0 && hkl.l() == 0 && hkl.k() == 0) { continue; } // Don't include the 0,0,0 reflection

        clipper::datatypes::F_phi<float> fbest_reflection = fbest[ih];
        clipper::datatypes::F_phi<float> fdiff_reflection = fdiff[ih];
        clipper::datatypes::F_sigF<float> fobs_reflection = m_fobs[ih];

        recalculated_data.emplace_back(hkl.h());
        recalculated_data.emplace_back(hkl.k());
        recalculated_data.emplace_back(hkl.l());
        recalculated_data.emplace_back(clipper::Util::rad2d(fobs_reflection.f()));
        recalculated_data.emplace_back(clipper::Util::rad2d(fobs_reflection.sigf()));
        recalculated_data.emplace_back(clipper::Util::rad2d(fbest_reflection.f()));
        recalculated_data.emplace_back(clipper::Util::rad2d(fbest_reflection.phi()));
        recalculated_data.emplace_back(clipper::Util::rad2d(fdiff_reflection.f()));
        recalculated_data.emplace_back(clipper::Util::rad2d(fdiff_reflection.phi()));
    }

    gemmi::Mtz new_mtz;
    new_mtz.set_cell_for_all(m_mtz.cell);
    new_mtz.set_spacegroup(m_mtz.spacegroup);
    new_mtz.add_dataset("SAILS");
    new_mtz.add_base();
    new_mtz.add_column("FP", 'F', -1, -1, true);
    new_mtz.add_column("SIGFP", 'Q', -1, -1, true);
    new_mtz.add_column("FWT", 'F', -1, -1, true);
    new_mtz.add_column("PHWT", 'P', -1, -1, true);
    new_mtz.add_column("DELFWT", 'F', -1, -1, true);
    new_mtz.add_column("PHDELWT", 'P', -1, -1, true);
    new_mtz.set_data(recalculated_data.data(), recalculated_data.size());
    new_mtz.ensure_asu();

    m_mtz = std::move(new_mtz);
    m_grid = load_grid(m_mtz, "FWT", "PHWT", false);
    m_difference_grid = load_grid(m_mtz, "DELFWT", "PHDELWT", true);
}

void Sails::XtalDensity::calculate_po_pc_map(gemmi::Structure &structure) {
    std::vector<clipper::Atom> atoms;
    form_atom_list(structure, atoms);

    clipper::Xmap<float> calculated_map = {m_spacegroup, m_cell, m_grid_sampling};
    clipper::EDcalc_iso<float> ed_calc = {m_resolution.limit()};
    clipper::HKL_data<clipper::data32::F_phi> fphic(m_hkl_info);

    ed_calc(calculated_map, atoms);
    calculated_map.fft_to(fphic);

    clipper::Xmap<float> new_map = m_best_map;
    new_map -= calculated_map;
    clipper::HKL_data<clipper::data32::F_phi> diff(m_fobs);
    new_map.fft_to(diff);

    std::vector<float> recalculated_data;
    int num_columns = 5;
    recalculated_data.reserve(m_fobs.data_size() * num_columns);
    for (HRI ih = m_fobs.first(); !ih.last(); ih.next()) {
        clipper::HKL hkl = ih.hkl();
        if (hkl.h() == 0 && hkl.l() == 0 && hkl.k() == 0) { continue; } // Don't include the 0,0,0 reflection

        clipper::datatypes::F_phi<float> diff_reflection = diff[ih];

        recalculated_data.emplace_back(hkl.h());
        recalculated_data.emplace_back(hkl.k());
        recalculated_data.emplace_back(hkl.l());
        recalculated_data.emplace_back(clipper::Util::rad2d(diff_reflection.f()));
        recalculated_data.emplace_back(clipper::Util::rad2d(diff_reflection.phi()));
    }

    gemmi::Mtz new_mtz;
    new_mtz.set_cell_for_all(m_mtz.cell);
    new_mtz.set_spacegroup(m_mtz.spacegroup);
    new_mtz.add_dataset("SAILS");
    new_mtz.add_base();
    new_mtz.add_column("FDIFFCALC", 'F', -1, -1, true);
    new_mtz.add_column("PDIFFCALC", 'P', -1, -1, true);
    new_mtz.set_data(recalculated_data.data(), recalculated_data.size());
    new_mtz.ensure_asu();

    m_po_pc_grid = load_grid(new_mtz, "FDIFFCALC", "PDIFFCALC", false);
}
