//
// Created by Jordan Dialpuri on 12/03/2025.
//

#ifndef SAILS_MATHS_H
#define SAILS_MATHS_H

#include <iostream>
#include <numeric>

namespace Sails::Maths {

    struct MeanAndVariance {
        MeanAndVariance() = default;
        MeanAndVariance(const double mean, const double variance): mean(mean), variance(variance) {}
        double mean;
        double variance;
    };

    MeanAndVariance calculate_mean_and_variance(const std::vector<double>& input);
}


#endif //SAILS_MATHS_H
