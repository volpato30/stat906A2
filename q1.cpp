//
// Created by rui on 28/10/18.
//

#include "q1.h"
#include "helper.h"

using std::default_random_engine;
using std::uniform_real_distribution;

default_random_engine generator;

double f_x(double x) {
    return 1.0 / (1 + pow(x, 2.0) + pow(x, 4.0));
}

double z_x(double x) {
    // Ez equals 28 / 27
    return 10 / 9.0 * pow(x, 2.0) - 32 / 9.0 * x + 28 / 9.0;
}

double optimalB(const std::vector<double> &x , const std::vector<double> &y) {
    unsigned long n = x.size();
    if (n != y.size()) {
        throw std::invalid_argument("x and y should have same length");
    }
    double mean_x = mean(x);
    double mean_y = mean(y);
    double numerator = 0;
    double denominator = 0;
    for (int i=0; i<n; i++) {
        numerator += (x[i] - mean_x) * (y[i] - mean_y);
        denominator += (x[i] - mean_x) * (x[i] - mean_x);
    }
    return numerator / denominator;
}


double vanillaMC1(unsigned long n) {
    // do a MC with n samples, return the estimated value.
    double mean = 0;
    uniform_real_distribution<double> rGen(0.0, 2.0);
    for (int i = 0; i < n; i++) {
        mean += f_x(rGen(generator));
    }
    return mean / n;
}

double controlVariateMC(unsigned long n, double b) {
    double output = 0;
    const double e_z = 28.0 / 27;
    uniform_real_distribution<double> rGen(0.0, 2.0);
    double temp;
    for (int i = 0; i < n; i++) {
        temp = rGen(generator);
        output += f_x(temp) - b * (z_x(temp) - e_z);
    }
    return output / n;
}


void q1() {
    // first do a pilot experiment to find optimal b
    unsigned long pilot_n = 10000;
    unsigned long number_of_experiments = 10000;
    unsigned long n = 100;
    std::vector<double> fx(pilot_n);
    std::vector<double> zx(pilot_n);
    uniform_real_distribution<double> rGen(0.0, 2.0);
    double temp;
    for (int i=0; i < int(pilot_n); i++) {
        temp = rGen(generator);
        fx[i] =f_x(temp);
        zx[i] = z_x(temp);
    }
    double b = optimalB(zx, fx);

    std::vector<double> vanillaMCOutput(number_of_experiments);
    std::vector<double> controlVariateMCOutput(number_of_experiments);
    for (int i=0; i < int(number_of_experiments); i++) {
        vanillaMCOutput[i] = vanillaMC1(n);
        controlVariateMCOutput[i] = controlVariateMC(n, b);
    }

    std::cout << "q1:" << std::endl;
    std::cout << "optimal b:\t" << b << std::endl;

    std::cout << "Vanilla MC:\t" << "mean:\t" << mean(vanillaMCOutput) << "\tvariance:\t" << computeVariance(vanillaMCOutput) << std::endl;
    std::cout << "control variate MC:\t" << "mean:\t" << mean(controlVariateMCOutput) << "\tvariance:\t" << computeVariance(controlVariateMCOutput) << std::endl;
    std::cout << "efficiency gain:\t" << computeVariance(vanillaMCOutput) / computeVariance(controlVariateMCOutput) << std::endl;
}