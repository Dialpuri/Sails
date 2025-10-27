//
// Created by Jordan Dialpuri on 22/10/2025.
//

#ifndef SAILS_SCORE_H
#define SAILS_SCORE_H
#include "sails-model.h"
#include "density/sails-density.h"

namespace Sails::Score {

    std::map<Glycosite, double> calculate_rsccs(Sails::Density* density, gemmi::Structure* structure, ResidueDatabase &residue_database);

    std::map<Glycosite, double> calculate_qscores(Sails::Density* density, gemmi::Structure* structure, ResidueDatabase &residue_database);

    double calculate_clash_score(Sails::Glycosite &site, gemmi::Structure* structure);

    namespace QScore {
        std::vector<gemmi::Position> fibonacci_sphere(int samples, float radius, const gemmi::Position &center);

        std::vector<gemmi::Position> get_radial_points(const gemmi::Position &position, float radius, int N, Glycosite& site, gemmi::NeighborSearch& ns);

        std::vector<double> sample_density(const gemmi::Grid<> *grid, std::vector<gemmi::Position>& positions);

        double calculate_q_score(const gemmi::Position & position, Glycosite &site, const gemmi::Grid<> *grid,
                                 gemmi::NeighborSearch &ns, float A, float B, float sigma, int N);
    }
}

#endif //SAILS_SCORE_H
