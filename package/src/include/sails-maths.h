//
// Created by Jordan Dialpuri on 12/03/2025.
//

#ifndef SAILS_MATHS_H
#define SAILS_MATHS_H

#include <iostream>
#include <numeric>
#include <gemmi/unitcell.hpp>

#include "sails-model.h"

namespace Sails::Maths {

    struct MeanAndVariance {
        MeanAndVariance() = default;
        MeanAndVariance(const double mean, const double variance): mean(mean), variance(variance) {}
        double mean;
        double variance;
    };
    using Matrix = std::vector<std::vector<double>>;

    MeanAndVariance calculate_mean_and_variance(const std::vector<float>& values);

    std::vector<gemmi::Position> fibonacci_sphere(int samples, float radius, const gemmi::Position& center);

    std::vector<double> linspace(double start, double end, int N);

    Matrix zeros(int rows, int cols);

    std::vector<double> fill(double value, int N);

    void assign_columns(Matrix& matrix, int col, const std::vector<double>& values);

    std::vector<double> row_mean(const Matrix& matrix);

    Matrix normalise_rows(const Matrix& matrix);

    double dot_product(const Matrix& A, const Matrix& B);

    double frobenius_norm(const Matrix& matrix);

    namespace Utils {
        std::vector<gemmi::Position> get_radial_points(gemmi::Position& pos, float radius, int N);
    }
}


#endif //SAILS_MATHS_H
