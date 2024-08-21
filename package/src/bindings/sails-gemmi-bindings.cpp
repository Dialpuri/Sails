//
// Created by Jordan Dialpuri on 22/07/2024.
//

#include "../include/sails-gemmi-bindings.h"
#include "gemmi/mtz.hpp"

#include <src/include/sails-utils.h>


gemmi::Mtz Sails::form_gemmi_mtz(MTZ& mtz) {
    gemmi::UnitCell cell = mtz.cell.to_gemmi_cell();
    gemmi::SpaceGroup spg = *gemmi::find_spacegroup_by_name(mtz.spacegroup);

    std::vector<float> data;
    int num_columns = 6;
    data.reserve(mtz.reflections.size() * num_columns);
    bool fwt_map_available = false;
    bool delfwt_map_available = false;

    for (const auto& reflection: mtz.reflections) {
        data.emplace_back(reflection.hkl.h);
        data.emplace_back(reflection.hkl.k);
        data.emplace_back(reflection.hkl.l);
        data.emplace_back(reflection.f_sigf.value1);
        data.emplace_back(reflection.f_sigf.value2);
        if (!std::isnan(reflection.fwt_phwt.value1)) {
            fwt_map_available = true;
            data.emplace_back(reflection.fwt_phwt.value1);
            data.emplace_back(reflection.fwt_phwt.value2);
        }
        if (!std::isnan(reflection.delfwt_phdelwt.value1)) {
            delfwt_map_available = true;
            data.emplace_back(reflection.delfwt_phdelwt.value1);
            data.emplace_back(reflection.delfwt_phdelwt.value2);
        }
    }
    gemmi::Mtz new_mtz = gemmi::Mtz(true);
    new_mtz.set_cell_for_all(cell);
    new_mtz.spacegroup_name = mtz.spacegroup;
    new_mtz.setup_spacegroup();
    new_mtz.add_dataset("SAILS");
    new_mtz.add_column("FP", 'F', -1, -1, true);
    new_mtz.add_column("SIGFP", 'Q', -1, -1, true);
    if (fwt_map_available) {
        new_mtz.add_column("FWT", 'F', -1, -1, true);
        new_mtz.add_column("PHWT", 'P', -1, -1, true);
    }
    if (delfwt_map_available) {
        new_mtz.add_column("DELFWT", 'F', -1, -1, true);
        new_mtz.add_column("PHDELWT", 'P', -1, -1, true);;
    }
    new_mtz.set_data(data.data(), data.size());
    new_mtz.update_reso();
    new_mtz.ensure_asu();
    return new_mtz;
}

Sails::MTZ Sails::form_sails_mtz(const gemmi::Mtz &mtz, const std::string &f_label, const std::string &sigf_label) {
    gemmi::Mtz::Column h = mtz.get_column_with_label("H");
    gemmi::Mtz::Column k = mtz.get_column_with_label("K");
    gemmi::Mtz::Column l = mtz.get_column_with_label("L");
    gemmi::Mtz::Column f = mtz.get_column_with_label(f_label);
    gemmi::Mtz::Column sigf = mtz.get_column_with_label(sigf_label);

    const gemmi::Mtz::Column* fwt = mtz.column_with_label("FWT");
    const gemmi::Mtz::Column* phwt = mtz.column_with_label("PHWT");

    const gemmi::Mtz::Column* delfwt = mtz.column_with_label("DELFWT");
    const gemmi::Mtz::Column* phdelwt = mtz.column_with_label("PHDELWT");

    if (h.size() != f.size()) {
        std::cout << h.size() << " " << k.size() << " " << l.size() << " " << f.size() << " " << sigf.size() << std::endl;
        throw std::runtime_error("Columns are unequal length");
    }

    std::vector<Reflection> reflections;
    reflections.reserve(h.size());
    for (int i = 0; i < h.size(); i++) {
        HKL hkl = {static_cast<int>(h[i]), static_cast<int>(k[i]), static_cast<int>(l[i])};
        NumberPair fsigf = {f[i], sigf[i]};

        NumberPair fwt_phwt;
        if (fwt != nullptr && phwt != nullptr) {
            fwt_phwt.value1 = fwt->operator[](i);
            fwt_phwt.value2 = phwt->operator[](i);
        }

        NumberPair delfwt_phdelwt;
        if (fwt != nullptr && phwt != nullptr) {
            delfwt_phdelwt.value1 = delfwt->operator[](i);
            delfwt_phdelwt.value2 = phdelwt->operator[](i);
        }

        Reflection reflection = {hkl, fsigf, fwt_phwt, delfwt_phdelwt};
        reflections.emplace_back(reflection);
    }

    Cell cell = Cell(mtz.cell);
    std::string spacegroup = mtz.spacegroup_name;
    return MTZ(reflections, cell, spacegroup);
}
