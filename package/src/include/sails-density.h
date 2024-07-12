//
// Created by Jordan Dialpuri on 07/07/2024.
//

#ifndef SAILS_SAILS_DENSITY_H
#define SAILS_SAILS_DENSITY_H

#include "sails-glycan.h"
#include <gemmi/mtz.hpp>
#include <gemmi/grid.hpp>
#include <gemmi/fourier.hpp>
#include <gemmi/dencalc.hpp>
#include "gemmi/it92.hpp"

namespace Sails {
	struct SuperpositionResult; // forward declaration from sails-linkage

	enum DensityScoreMethod {
		atomwise, rscc, rsr
	};


	class Density {
	public:
		explicit Density(const std::string& mtz_path, const std::string& f_col, const std::string& phi_col);
		explicit Density(gemmi::Mtz& mtz, const std::string& f_col, const std::string& phi_col);

		double score_residue(gemmi::Residue &residue, const DensityScoreMethod &method = atomwise);

	// private:

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
		gemmi::Grid<> load_grid(const gemmi::Mtz &mtz, const std::string& f_col,
		                                            const std::string& phi_col, bool normalise);


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
		[[nodiscard]] float atomwise_score(const gemmi::Residue& residue) const;

		gemmi::Grid<> calculate_density_for_box(gemmi::Residue &residue);

		float rscc_score(gemmi::Residue &residue);

		float calculate_rscc(std::vector<float> obs_values, std::vector<float> calc_values);

		float rscc_score(SuperpositionResult& result);

		float rsr_score(gemmi::Residue &residue);

		void initialise_density_calculator() {
			density_calculator.grid.copy_metadata_from(m_grid);
			density_calculator.d_min = m_mtz.resolution_high();
			density_calculator.initialize_grid();
		}


	private:
		gemmi::DensityCalculator<gemmi::IT92<float>, float > density_calculator;
		gemmi::Grid<> m_grid{};
		gemmi::Mtz m_mtz;

		std::map<std::string, gemmi::Grid<>> calculated_maps;
	};

}// namespace Sails

#endif //SAILS_SAILS_DENSITY_H
