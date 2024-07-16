//
// Created by Jordan Dialpuri on 07/07/2024.
//

#include "../include/sails-density.h"

#include <clipper/contrib/edcalc.h>
#include <clipper/contrib/sfweight.h>
#include <clipper/minimol/minimol.h>

#include "../include/sails-refine.h"


Sails::Density::Density(const std::string &mtz_path) {
    gemmi::Mtz mtz = gemmi::read_mtz_file(mtz_path);
    m_mtz = std::move(mtz);
    initialise_hkl();
}

Sails::Density::Density(gemmi::Mtz &mtz) {
    m_mtz = std::move(mtz);
    initialise_hkl();
}

double Sails::Density::score_residue(gemmi::Residue &residue, const DensityScoreMethod &method) {
    switch (method) {
        case atomwise:
            return atomwise_score(residue);
        case rscc:
            return rscc_score(residue);
        case rsr:
            return rsr_score(residue);
        default:
            return -1;
    }
}

gemmi::Grid<> Sails::Density::load_grid(const gemmi::Mtz &mtz, const std::string &f_col, const std::string &phi_col,
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

void Sails::Density::initialise_hkl() {
    m_resolution = clipper::Resolution(m_mtz.resolution_high());
    const gemmi::UnitCell *cell = &m_mtz.cell;
    m_cell = clipper::Cell(clipper::Cell_descr(
        cell->a, cell->b, cell->c, cell->alpha, cell->beta, cell->gamma
    ));
    m_spacegroup = clipper::Spacegroup(clipper::Spgr_descr(m_mtz.spacegroup_name));
    m_hkl_info = {m_spacegroup, m_cell, m_resolution, true};
    m_grid_sampling = {m_spacegroup, m_cell, m_resolution};
}

void Sails::Density::load_hkl(const std::string &f, const std::string &sig_f) {
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

void Sails::Density::form_atom_list(const gemmi::Structure &structure, std::vector<clipper::Atom>& atoms) {
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

void Sails::Density::recalculate_map(gemmi::Structure &structure) {
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

    std::vector<float> recalculated_data;
    int num_columns = 6;
    recalculated_data.reserve(fdiff.data_size() * num_columns);
    for (HRI ih = m_fobs.first(); !ih.last(); ih.next() ) {
        clipper::HKL hkl = ih.hkl();
        if (hkl.h() == 0 && hkl.l() == 0 && hkl.k() == 0) { continue;} // Don't include the 0,0,0 reflection

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
    new_mtz.add_column("PHDELWT", 'P', -1, -1, true);;
    new_mtz.set_data(recalculated_data.data(), recalculated_data.size());
    new_mtz.ensure_asu();

    m_mtz = std::move(new_mtz);
    m_grid = load_grid(m_mtz, "FWT", "PHWT", false);
    m_difference_grid = load_grid(m_mtz, "DELFWT", "PHDELWT", true);
}

void Sails::Density::write_mtz(const std::string &path) const {
    m_mtz.write_to_file(path);
}

float Sails::Density::atomwise_score(const gemmi::Residue &residue) const {
    return std::transform_reduce(residue.atoms.begin(), residue.atoms.end(), 0.0f, std::plus<>(),
                                 [&](const gemmi::Atom &current_atom) {
                                     return m_grid.interpolate_value(current_atom.pos);
                                 });
}

gemmi::Grid<> Sails::Density::calculate_density_for_box(gemmi::Residue &residue) const {
    gemmi::DensityCalculator<gemmi::IT92<float>, float> density_calculator;
    density_calculator.grid.copy_metadata_from(m_grid);
    density_calculator.d_min = m_mtz.resolution_high();
    density_calculator.initialize_grid();
    for (auto &atom: residue.atoms) {
        density_calculator.add_atom_density_to_grid(atom);
    }
    density_calculator.grid.symmetrize_sum();
    return density_calculator.grid;
}

float Sails::Density::calculate_rscc(std::vector<float> obs_values, std::vector<float> calc_values) {
    if (obs_values.size() != calc_values.size())
        throw std::runtime_error(
            "RSCC obs and calc lists are different sizes");

    if (obs_values.empty()) throw std::runtime_error("Observation list is empty");
    if (calc_values.empty()) throw std::runtime_error("Calculated list is empty");

    float obs_average = std::accumulate(obs_values.begin(), obs_values.end(), 0.0f) / obs_values.size();
    float calc_average = std::accumulate(calc_values.begin(), calc_values.end(), 0.0f) / calc_values.size();

    if (calc_average == 0.0f) throw std::runtime_error("Calculated map average is 0");

    float numerator = 0.0f;
    float obs_sum_sq = 0.0f;
    float calc_sum_sq = 0.0f;

    for (int i = 0; i < obs_values.size(); i++) {
        float obs_delta = obs_values[i] - obs_average;
        float calc_delta = calc_values[i] - calc_average;

        numerator += obs_delta * calc_delta;
        obs_sum_sq += (obs_delta * obs_delta);
        calc_sum_sq += (calc_delta * calc_delta);
    }

    float denominator = sqrt(obs_sum_sq * calc_sum_sq);

    if (denominator == 0.0f) throw std::runtime_error("RSCC Denominator is 0");
    return numerator / denominator;
}


float Sails::Density::rscc_score(gemmi::Residue &residue) {
    if (residue.atoms.empty()) throw std::runtime_error("Residue is empty during RSCC check");

    gemmi::Box<gemmi::Position> box;
    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }
    box.add_margin(1);

    gemmi::Grid<> calc = calculate_density_for_box(residue);

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    std::vector<float> obs_values = {};
    std::vector<float> calc_values = {};

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                obs_values.emplace_back(m_grid.interpolate_value(position));
                calc_values.emplace_back(calc.interpolate_value(position));
            }
        }
    }

    return calculate_rscc(obs_values, calc_values);
}

float Sails::Density::rscc_score(SuperpositionResult &result) {
    gemmi::Box<gemmi::Position> box;
    gemmi::Residue residue = result.new_residue;

    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }
    box.add_margin(1);

    if (calculated_maps.find(residue.name) == calculated_maps.end()) {
        gemmi::Grid<> reference = calculate_density_for_box(result.reference_residue);
        calculated_maps[residue.name] = std::move(reference);
    }
    gemmi::Grid<> *calculated = &calculated_maps[residue.name];

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    std::vector<float> obs_values = {};
    std::vector<float> calc_values = {};

    gemmi::Residue r1, r2, r3;

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                obs_values.emplace_back(m_grid.interpolate_value(position));
                gemmi::Vec3 translated_position = result.transformation.inverse().apply(position);
                calc_values.emplace_back(calculated->interpolate_value(gemmi::Position(translated_position)));
            }
        }
    }

    return calculate_rscc(obs_values, calc_values);
}

float Sails::Density::rsr_score(gemmi::Residue &residue) {
    gemmi::Box<gemmi::Position> box;
    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }

    gemmi::Grid<> calc_grid = calculate_density_for_box(residue);

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    float numerator = 0.0f;
    float denominator = 0.0f;

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                float obs = m_grid.interpolate_value(position);
                float calc = calc_grid.interpolate_value(position);
                numerator += abs(obs - calc);
                denominator += abs(obs + calc);
            }
        }
    }

    if (denominator == 0.0f) throw std::runtime_error("Box is empty");
    return numerator / denominator;
}

float Sails::Density::rsr_score(SuperpositionResult &result) {
    gemmi::Box<gemmi::Position> box;
    gemmi::Residue residue = result.new_residue;

    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }
    box.add_margin(1);

    // calculate map if not found
    if (calculated_maps.find(residue.name) == calculated_maps.end()) {
        gemmi::Grid<> reference = calculate_density_for_box(result.reference_residue);
        calculated_maps[residue.name] = std::move(reference);
    }
    gemmi::Grid<> *calculated = &calculated_maps[residue.name];

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    float numerator = 0.0f;
    float denominator = 0.0f;

    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                float obs = m_grid.interpolate_value(position);
                gemmi::Vec3 translated_position = result.transformation.inverse().apply(position);
                float calc = calculated->interpolate_value(gemmi::Position(translated_position));
                numerator += abs(obs - calc);
                denominator += abs(obs + calc);
            }
        }
    }
    if (denominator == 0.0f) throw std::runtime_error("Box is empty");
    return numerator / denominator;
}

float Sails::Density::difference_density_score(gemmi::Residue &residue) const {
    gemmi::Box<gemmi::Position> box;
    for (auto &atom: residue.atoms) {
        box.extend(atom.pos);
    }

    const gemmi::Position max = box.maximum;
    const gemmi::Position min = box.minimum;

    float sum = 0.0f;
    int points = 0;
    constexpr double step_size = 0.5;
    for (double x = min.x; x <= max.x; x += step_size) {
        for (double y = min.y; y <= max.y; y += step_size) {
            for (double z = min.z; z <= max.z; z += step_size) {
                gemmi::Position position = {x, y, z};
                float value = m_difference_grid.interpolate_value(position);
                sum += abs(value);
                points++;
            }
        }
    }

    return sum/points;
}

float Sails::Density::score_atomic_position(const gemmi::Atom &atom) const {
    return score_position(atom.pos);
}


float Sails::Density::score_position(const gemmi::Position &pos) const {
    return m_grid.interpolate_value(pos);
}

