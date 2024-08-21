//
// Created by Jordan Dialpuri on 13/08/2024.
//

#include "sails-density.h"
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>

namespace Sails {

    class XtalDensity : public Density {
    public:
        explicit XtalDensity(gemmi::Mtz &mtz);

        [[nodiscard]] const gemmi::Mtz *get_mtz() const override { return &m_mtz; }

        [[nodiscard]] const gemmi::Grid<> *get_work_grid() const override { return &m_po_pc_grid; }

        [[nodiscard]] const gemmi::Grid<> *get_best_grid() const override { return &m_grid; }

        [[nodiscard]] const gemmi::Grid<> *get_difference_grid() const override { return &m_difference_grid; }

        [[nodiscard]] const double get_resolution() const override { return 2.0; }

        [[nodiscard]] const DensityScoreMethod get_score_method() const override { return score_method;}

        [[nodiscard]] std::unordered_map <std::string, gemmi::Grid<>> *get_calculated_maps() override {
            return &calculated_maps;
        }

        /**
         * @brief Recalculates the map based on the given structure.
         *
         * This method takes a gemmi::Structure object and recalculates the map with a clipper sigma-a calculation.
         * This function updates the internal mtz, grid and density grid.
         *
         * @param structure The gemmi::Structure object for which the map needs to be recalculated.
         */
        void recalculate_map(gemmi::Structure &structure);

        /**
		* @brief Recalculates the map based on the given structure.
		*
		* This method takes a gemmi::Structure object and recalculates the p_o - p_c map, where p_o is the
		* internal best map and p_c is calculated with the provided structure. This function should be called after
		* recalculate_map
		* This function updates the internal mtz, grid and density grid.
		*
		* @param structure The gemmi::Structure object for which the map needs to be recalculated.
		*/
		void calculate_po_pc_map(gemmi::Structure &structure);

    private:
        /**
		 * @brief Initializes the necessary parameters for density calculation.
		 *
		 * This method initializes the necessary parameters for density calculation, such as resolution, cell, spacegroup, hkl_info, and grid_sampling.
		 * These parameters are obtained from the m_mtz object which represents the input MTZ file.
		 */
        void initialise_hkl();

        /**
         * @brief Loads the reflection data from the specified MTZ file and builds HKL data.
         *
         * This method reads the reflection data from the MTZ file specified by 'f' and 'sig_f' and builds
         * HKL (H, K, L) data. It populates the internal structures with the HKL data for further analysis.
         *
         * @param f The path to the MTZ file containing the reflection data.
         * @param sig_f The path to the MTZ file containing the sigma data.
         */
        void load_hkl(const std::string &f, const std::string &sig_f);

        /**
		 * @brief Transforms F and Phi values from an MTZ file into a grid map.
		 *
		 * This method takes F and Phi values from an MTZ file and transforms them into a grid using an FFT.
		 *
		 * @param mtz The MTZ file containing the F and Phi values.
		 * @param f_col The label of the column in the MTZ file containing the F values e.g. FWT.
		 * @param phi_col The label of the column in the MTZ file containing the Phi values e.g PHWT.
		 * @param normalise Indicates whether to normalize the resulting grid map.
		 * @return The real space grid formed after FFT.
		 */
        static gemmi::Grid<> load_grid(const gemmi::Mtz &mtz, const std::string &f_col,
                                       const std::string &phi_col, bool normalise);

        /**
         * @brief Forms a list of clipper::Atom objects from a gemmi::Structure object.
         *
         * This method takes a gemmi::Structure object and populates a vector of clipper::Atom objects.
         * Each atom in the structure is converted to a clipper::Atom object and added to the list.
         *
         * @param structure The gemmi::Structure object from which to form the list of atoms.
         * @param atoms The vector of clipper::Atom objects to be populated.
         */
        static void form_atom_list(const gemmi::Structure &structure, std::vector <clipper::Atom> &atoms);

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
         * @brief Density score method.
         *
         * The DensityScoreMethod class is used to represent the score method for scoring residues to density
         */
        DensityScoreMethod score_method = atomwise;

        /**
         * Fc maps for residues in standard positions - used for fast RSCC calculations
         */
        std::unordered_map <std::string, gemmi::Grid<>> calculated_maps;

        // clipper HKL - initialised in initialise_hkl
        /**
         * Clipper Spacegroup used for map recalculation - initialised in initialise_hkl
         */
        clipper::Spacegroup m_spacegroup;

        /**
         * Clipper Resolution used for map recalculation - initialised in initialise_hkl
         */
        clipper::Resolution m_resolution;

        /**
         * Clipper Cell used for map recalculation - initialised in initialise_hkl
         */
        clipper::Cell m_cell;

        /**
         * Clipper container for reflection metadata used for map recalculation
         */
        clipper::HKL_info m_hkl_info;

        /**
         * Clipper Grid sampling used for map recalcualtion - initialised in initialise_hkl
         */
        clipper::Grid_sampling m_grid_sampling;

        /**
         * Clipper list of F obs reflections
         */
        clipper::HKL_data <clipper::data32::F_sigF> m_fobs; // initialsised in load_hkl

        /**
         * Clipper best map
         */
        clipper::Xmap<float> m_best_map;
    };
}