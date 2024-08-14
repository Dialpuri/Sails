//
// Created by Jordan Dialpuri on 07/07/2024.
//

#ifndef SAILS_SAILS_DENSITY_H
#define SAILS_SAILS_DENSITY_H



#include <clipper/clipper.h>
#include <gemmi/mtz.hpp>
#include <gemmi/grid.hpp>
#include <gemmi/dencalc.hpp>
#include <gemmi/fourier.hpp>
#include <gemmi/it92.hpp>
#include <gemmi/calculate.hpp>
#include <gemmi/ccp4.hpp>
#include <gemmi/modify.hpp>


namespace Sails {
	struct SuperpositionResult;
	typedef clipper::HKL_info::HKL_reference_index HRI;

	enum DensityScoreMethod {
		atomwise, rscc, rsr, dds
	};

	class Density {
    public:
        Density() = default;

        [[nodiscard]] virtual const gemmi::Mtz* get_mtz() const = 0;

        [[nodiscard]] virtual const gemmi::Grid<>* get_work_grid() const = 0;

        [[nodiscard]] virtual const gemmi::Grid<>* get_best_grid() const = 0;

        [[nodiscard]] virtual const gemmi::Grid<>* get_difference_grid() const = 0;

        [[nodiscard]] virtual const double get_resolution() const = 0;

        [[nodiscard]] virtual std::unordered_map<std::string, gemmi::Grid<>>* get_calculated_maps()  = 0;


        /**
         * @brief Scores a residue based on the specified density score method.
         *
         * This method takes a gemmi::Residue object and scores it based on the specified density score method.
         *
         * @param residue The gemmi::Residue object to be scored.
         * @param method The density score method to be used. Default value is atomwise.
         *
         * @return The score of the residue based on the specified density score method.
         */
        double score_residue(gemmi::Residue &residue, const DensityScoreMethod &method = atomwise);

        /**
		 * @brief Calculates the atomwise score for a given residue.
		 *
		 * This method calculates the atomwise score by iteratively interpolating the grid values at each atom position
		 * in the given residue and summing them up.
		 *
		 * @param residue The residue for which the atomwise score is calculated.
		 *
		 * @return The atomwise score for the given residue.
		 */
        [[nodiscard]] float atomwise_score(const gemmi::Residue &residue) const;

        /**
         * @brief Calculates the density for a given box based on a gemmi::Residue object.
         *
         * This method takes a gemmi::Residue object and calculates the density for the specified box
         * using the gemmi::DensityCalculator class. The density calculation is performed using the
         * density score method specified in the constructor of the gemmi::DensityCalculator.
         *
         * @param residue The gemmi::Residue object for which the density is calculated.
         *
         * @return The calculated density grid for the specified box.
         */
        gemmi::Grid<> calculate_density_for_box(gemmi::Residue &residue) const;

        /**
         * @brief Calculates the RSCC (Real Space Correlation Coefficient) score for a given residue.
         *
         * The RSCC score is a measure of how well the observed electron density matches the calculated electron density
         * for a residue. This method calculates the RSCC score by comparing the observed and calculated density values
         * at different positions within a bounding box around the residue.
         *
         * @param result The gemmi::Residue object for which the RSCC score will be calculated.
         *
         * @return The RSCC score for the given residue.
         * @throws std::runtime_error if the residue is empty.
         */
        float rscc_score(gemmi::Residue& residue) const;

        /**
         * @brief Calculates the Real Space Correlation Coefficient (RSCC) between observed and calculated values.
         *
         * This method takes in two vectors of observed and calculated values and calculates the RSCC.
         *
         * @param obs_values The vector of observed values.
         * @param calc_values The vector of calculated values.
         * @throws std::runtime_error If the sizes of the observation and calculation lists are different, if either list is empty, if the calculated map average is 0, or if the denominator is 0.
         *
         * @return The RSCC between the observed and calculated values.
         */
        static float calculate_rscc(std::vector<float> obs_values, std::vector<float> calc_values) ;

        /**
         * @brief Calculates the RSCC score for a given superposition result.
         *
         * This method takes a SuperpositionResult object and calculates the RSCC (Real Space Correlation Coefficient) score
         * for it based on the following steps:
         * - Calculates the minimum bounding box for the new_residue atoms.
         * - Extends the box by adding a margin of 1 Angstrom.
         * - Checks if the density map for the residue name is already calculated. If not, it calculates it using the
         *   calculate_density_for_box() method and stores it in the calculated_maps.
         * - Retrieves the calculated density map for the residue name from the calculated_maps.
         * - Iterates over a grid with a step size of 0.5 (Angstrom or any other unit) within the bounding box and
         *   retrieves the observed and calculated density values at each position.
         * - Calculates and returns the RSCC score using the calculate_rscc() method with the observed and calculated density
         *   values.
         *
         * @param result The SuperpositionResult object for which the RSCC score needs to be calculated.
         *
         * @return The RSCC score for the given superposition result.
         */
        float rscc_score(SuperpositionResult &result);

        /**
         * @brief Calculates the RSR (Real Space R) score for a given residue.
         *
         * The RSR score measures the agreement between the observed electron density and the calculated electron density for a residue.
         * The score is calculated by comparing the observed and calculated electron densities at regular intervals within the residue's bounding box.
         * The RSR score is defined as the absolute difference between the observed and calculated densities divided by
         * the sum of the absolute values of the observed and calculated densities.
         *
         * @param residue The gemmi::Residue object for which the RSR score is calculated.
         *
         * @return The RSR score for the given residue.
         *
         * @throw std::runtime_error if the residue's bounding box is empty.
         */
        float rsr_score(gemmi::Residue &residue);

        /**
         * @brief Calculates the RSR (Real Space R) score for a given SuperpositionResult.
         *
         * This method calculates the RSR score for a given SuperpositionResult.
         * The RSR score measures the agreement between the observed electron density and the calculated electron density for a residue.
         * The score is calculated by comparing the observed and calculated electron densities at regular intervals within the residue's bounding box.
         * The RSR score is defined as the absolute difference between the observed and calculated densities divided by
         * the sum of the absolute values of the observed and calculated densities.
         *
         * @param result The SuperpositionResult object containing the required data for calculating the RSR score.
         *
         * @return The RSR score for the given SuperpositionResult.
         * @throws std::runtime_error if the box is empty (denominator is 0.0f).
         */
        float rsr_score(SuperpositionResult &result);

        /**
         * @brief Calculates the difference density score for a residue.
         *
         * This method calculates the difference density score for the given residue using the difference_grid.
         *
         * @param residue The gemmi::Residue object for which the difference density score is to be calculated.
         *
         * @return The difference density score for the residue.
         */
        float difference_density_score(gemmi::Residue &residue) const;

        /**
         * @brief Scores an atom
         *
         * This method takes a gemmi::Atom object and generates a score corresponding to the interpolated density value.
         *
         * @param atom The gemmi::Atom object whose position needs to be scored.
         *
         * @return The score of the atom based on the density value at that atom's position.
         */
        [[nodiscard]] float score_atomic_position(const gemmi::Atom& atom) const;

        /**
         * @brief Scores a position based on the density value at that position.
         *
         * This method takes a gemmi::Position object and scores it based on the density value at that position.
         *
         * @param pos The gemmi::Position object to be scored.
         *
         * @return The score of the position based on the density value at that position.
         */
        [[nodiscard]] float score_position(const gemmi::Position& pos) const;

    };

//    class Density {
//	public:
//		Density() = default;
//
//		explicit Density(const std::string &mtz_path);
//
//		explicit Density(gemmi::Mtz &mtz);
//
//		Density(Density &d) {
//			this->m_mtz = std::move(d.m_mtz);
//			this->m_grid = d.m_grid;
//			this->calculated_maps = d.calculated_maps;
//		};
//
//
//		// private:
//
//		/**
//		 * @brief Transforms F and Phi values from an MTZ file into a grid map.
//		 *
//		 * This method takes F and Phi values from an MTZ file and transforms them into a grid using an FFT.
//		 *
//		 * @param mtz The MTZ file containing the F and Phi values.
//		 * @param f_col The label of the column in the MTZ file containing the F values e.g. FWT.
//		 * @param phi_col The label of the column in the MTZ file containing the Phi values e.g PHWT.
//		 * @param normalise Indicates whether to normalize the resulting grid map.
//		 * @return The real space grid formed after FFT.
//		 */
//		static gemmi::Grid<> load_grid(const gemmi::Mtz &mtz, const std::string &f_col,
//		                               const std::string &phi_col, bool normalise);
//
//
//		/**
//		 * @brief Initializes the necessary parameters for density calculation.
//		 *
//		 * This method initializes the necessary parameters for density calculation, such as resolution, cell, spacegroup, hkl_info, and grid_sampling.
//		 * These parameters are obtained from the m_mtz object which represents the input MTZ file.
//		 */
//		void initialise_hkl();
//
//		/**
//		 * @brief Loads the reflection data from the specified MTZ file and builds HKL data.
//		 *
//		 * This method reads the reflection data from the MTZ file specified by 'f' and 'sig_f' and builds
//		 * HKL (H, K, L) data. It populates the internal structures with the HKL data for further analysis.
//		 *
//		 * @param f The path to the MTZ file containing the reflection data.
//		 * @param sig_f The path to the MTZ file containing the sigma data.
//		 */
//		void load_hkl(const std::string &f, const std::string &sig_f);
//
//		/**
//		 * @brief Forms a list of clipper::Atom objects from a gemmi::Structure object.
//		 *
//		 * This method takes a gemmi::Structure object and populates a vector of clipper::Atom objects.
//		 * Each atom in the structure is converted to a clipper::Atom object and added to the list.
//		 *
//		 * @param structure The gemmi::Structure object from which to form the list of atoms.
//		 * @param atoms The vector of clipper::Atom objects to be populated.
//		 */
//		static void form_atom_list(const gemmi::Structure &structure, std::vector<clipper::Atom> &atoms);
//
//		/**
//		 * @brief Recalculates the map based on the given structure.
//		 *
//		 * This method takes a gemmi::Structure object and recalculates the map with a clipper sigma-a calculation.
//		 * This function updates the internal mtz, grid and density grid.
//		 *
//		 * @param structure The gemmi::Structure object for which the map needs to be recalculated.
//		 */
//		void recalculate_map(gemmi::Structure &structure);
//
//		/**
//		* @brief Recalculates the map based on the given structure.
//		*
//		* This method takes a gemmi::Structure object and recalculates the p_o - p_c map, where p_o is the
//		* internal best map and p_c is calculated with the provided structure. This function should be called after
//		* recalculate_map
//		* This function updates the internal mtz, grid and density grid.
//		*
//		* @param structure The gemmi::Structure object for which the map needs to be recalculated.
//		*/
//		void calculate_po_pc_map(gemmi::Structure &structure);
//
//		/**
//		 * @brief Writes the density data to an MTZ file.
//		 *
//		 * This method writes the density data to an MTZ file at the specified path.
//		 *
//		 * @param path The path of the MTZ file to be written.
//		 */
//		void write_mtz(const std::string &path) const;
//
//
//
//		// private:
//		/**
//		 * Best map
//		 */
//		gemmi::Grid<> m_grid{};
//
//		/**
//		 * Difference map
//		 */
//		gemmi::Grid<> m_difference_grid{};
//
//		/**
//		 * Po-Pc map
//		 */
//		gemmi::Grid<> m_po_pc_grid{};
//
//		/**
//		 * MTZ Object
//		 */
//		gemmi::Mtz m_mtz;
//
//		/**
//		 * Fc maps for residues in standard positions - used for fast RSCC calculations
//		 */
//		std::unordered_map<std::string, gemmi::Grid<> > calculated_maps;
//
//		// clipper HKL - initialised in initialise_hkl
//		/**
//		 * Clipper Spacegroup used for map recalcualtion - initialised in initialise_hkl
//		 */
//		clipper::Spacegroup m_spacegroup;
//
//		/**
//		 * Clipper Resolution used for map recalcualtion - initialised in initialise_hkl
//		 */
//		clipper::Resolution m_resolution;
//
//		/**
//		 * Clipper Cell used for map recalcualtion - initialised in initialise_hkl
//		 */
//		clipper::Cell m_cell;
//
//		/**
//		 * Clipper container for reflection metadata used for map recalcualtion
//		 */
//		clipper::HKL_info m_hkl_info;
//
//		/**
//		 * Clipper Grid sampling used for map recalcualtion - initialised in initialise_hkl
//		 */
//		clipper::Grid_sampling m_grid_sampling;
//
//		/**
//		 * Clipper list of F obs reflections
//		 */
//		clipper::HKL_data<clipper::data32::F_sigF> m_fobs; // initialsised in load_hkl
//
//		/**
//		 * Clipper best map
//		 */
//		clipper::Xmap<float> m_best_map;
//	};
} // namespace Sails

#endif //SAILS_SAILS_DENSITY_H
