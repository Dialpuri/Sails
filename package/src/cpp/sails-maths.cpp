//
// Created by Jordan Dialpuri on 12/03/2025.
//


#include "../include/sails-maths.h"



Sails::Maths::MeanAndVariance Sails::Maths::calculate_mean_and_variance(const std::vector<double>& values) {
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
