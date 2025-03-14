//
// Created by Jordan Dialpuri on 12/03/2025.
//


#include "../include/sails-maths.h"

#include <clipper/core/clipper_types.h>

Sails::Maths::MeanAndVariance Sails::Maths::calculate_mean_and_variance(const std::vector<float>& values) {
    if (values.empty()) {
        return {0.0, 0.0};
    }

    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    double mean = sum / values.size();

    double variance = 0.0;
    for (const double value : values) {
        variance += (value - mean) * (value - mean);
    }
    variance /= values.size();

    return {mean, variance};
}


std::vector<gemmi::Position> Sails::Maths::fibonacci_sphere(const int samples, const float radius, const gemmi::Position& center) {
    std::vector<gemmi::Position> points;

    const float offset = 2.0 / samples;
    const float increment = M_PI * (3.0 - sqrt(5.0));
    for (int i = 0; i < samples; i++) {
        const float y = ((i * offset) - 1) + (offset / 2);
        const float r = sqrt(1 - y * y);

        const float phi = i * increment;

        const float x = cos(phi) * r;
        const float z = sin(phi) * r;

        gemmi::Position point = {x, y, z};
        point *= radius;
        point += center;
        points.emplace_back(point);
    }
    return points;
}

std::vector<double> Sails::Maths::linspace(const double start, const double end, const int N) {
    std::vector<double> result(N);
    const double step = (end - start) / (N - 1);

    for (int i = 0; i < N; ++i) {
        result[i] = start + i * step;
    }

    return result;
}

std::vector<std::vector<double>> Sails::Maths::zeros(const int rows, const int cols) {
    return std::vector<std::vector<double>>(rows, std::vector<double>(cols, 0.0));
}

std::vector<double> Sails::Maths::fill(const double value, const int N) {
    return std::vector<double>(N, value);
}

void Sails::Maths::assign_columns(std::vector<std::vector<double>> &matrix, const int col, const std::vector<double> &values) {
    const int rows = matrix.size();
    for (int i = 0; i < rows; ++i) {
        matrix[i][col] = values[i];
    }
}

std::vector<double> Sails::Maths::row_mean(const Matrix &matrix) {
    std::vector<double> means(matrix.size(), 0.0);

    for (size_t i = 0; i < matrix.size(); ++i) {
        means[i] = std::accumulate(matrix[i].begin(), matrix[i].end(), 0.0) / matrix[i].size();
    }

    return means;
}


Sails::Maths::Matrix Sails::Maths::normalise_rows(const Matrix &matrix) {
    std::vector<std::vector<double>> normalized = matrix;
    std::vector<double> means = row_mean(matrix);

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            normalized[i][j] -= means[i];
        }
    }

    return normalized;
}

double Sails::Maths::dot_product(const Matrix &A, const Matrix &B) {
    double sum = 0.0;

    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            sum += A[i][j] * B[i][j];
        }
    }

    return sum;
}

double Sails::Maths::frobenius_norm(const Matrix &matrix) {
    double sum = 0.0;
    for (const auto& row : matrix) {
        for (const double val : row) {
            sum += val * val;
        }
    }

    return std::sqrt(sum);
}


std::vector<gemmi::Position> Sails::Maths::Utils::get_radial_points(gemmi::Position &pos, float radius, int N) {
    std::vector<gemmi::Position> points;

    constexpr int max_iter = 200;
    for (int i = 0; i < max_iter; i++) {
        auto sampled_sphere = fibonacci_sphere(N, radius, pos);

        for (auto point : sampled_sphere) {
            points.emplace_back(point);
            if (points.size() >= N) {break;}
        }
        if (points.size() >= N) {break;}
    }

    return points;
}
