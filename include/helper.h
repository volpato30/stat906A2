//
// Created by rui on 28/10/18.
//

#ifndef STAT906A2_HELPER_H
#define STAT906A2_HELPER_H

#include <vector>
#include <cmath>

inline double mean(const std::vector<double> &x) {
    double m = 0;
    auto length = double(x.size());
    for (const double &d : x) {
        m += d;
    }
    m = m / length;
    return m;
}

inline double computeVariance(const std::vector<double> &x_array) {
    double mean = 0, variance = 0;
    unsigned long n = x_array.size();
    for (const double &x : x_array) {
        mean += x;
    }
    mean = mean / n;
    for (const double &x : x_array) {
        variance += pow(x - mean, 2);
    }
    variance = variance / (n - 1);
    return variance;
}

#endif //STAT906A2_HELPER_H
