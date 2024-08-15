//
// Created by Jordan Dialpuri on 13/08/2024.
//

#include "sails-density.h"

namespace Sails {

    class EMDensity : public Density {
    public:
        explicit EMDensity(gemmi::Grid<>& grid);

        [[nodiscard]] const gemmi::Mtz *get_mtz() const override { return &m_mtz; }

        [[nodiscard]] const gemmi::Grid<> *get_work_grid() const override { return &m_grid; }

        [[nodiscard]] const gemmi::Grid<> *get_best_grid() const override { return &m_grid; }

        [[nodiscard]] const gemmi::Grid<> *get_difference_grid() const override { return &m_grid; }

        [[nodiscard]] const double get_resolution() const override { return 2.0; }

        [[nodiscard]] std::unordered_map <std::string, gemmi::Grid<>> *get_calculated_maps() override {
            return &calculated_maps;
        }

    private:
        /**
		 * Best map
		 */
        gemmi::Grid<> m_grid{};

        /**
         * Difference map
         */
        gemmi::Grid<> m_difference_grid{};

        /**
         * Po-Pc map
         */
        gemmi::Grid<> m_po_pc_grid{};

        /**
         * MTZ Object
         */
        gemmi::Mtz m_mtz;

        /**
         * Fc maps for residues in standard positions - used for fast RSCC calculations
         */
        std::unordered_map <std::string, gemmi::Grid<>> calculated_maps{};

    };
}