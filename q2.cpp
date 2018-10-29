//
// Created by rui on 28/10/18.
//
#include "helper.h"
#include "q2.h"

using std::uniform_real_distribution;
using std::normal_distribution;

std::default_random_engine generator2;

double NormalCDFInverse(double p);
double RationalApproximation(double t);

double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) /
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p
           << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument( os.str() );
    }

    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}

void randomNormal(std::vector<double> &x) {
    std::normal_distribution<double> standard_normal{0,1};
    for (double &v : x) {
        v = standard_normal(generator2);
    }
}


double g_x(const std::vector<double> &x) {
    double result = 0;
    result = 2 * x[0] + 2 * x[1] * x[1] + x[2];
    result = sin(result);
    return result;
}

double vanillaMC2(unsigned long n) {
    // do a MC with n samples, return the estimated value.
    std::vector<double> z_vector(3);
    double mean = 0;
    for (int i = 0; i < n; i++) {
        randomNormal(z_vector);
        mean += g_x(z_vector);
    }
    return mean / n;
}

double kFoldStratify(unsigned int k, unsigned long n) {
    double mean = 0;
    auto num_rep = int(n / k);
    std::vector<double> x_vector(3);
    std::vector<double> z_vector(3);
    double y;
    uniform_real_distribution<double> runif{0.0, 1.0};
    for (int i = 1; i < k + 1; i++) {
        for (int j = 0; j < num_rep; j++) {
            y = (i - 1 + runif(generator2)) / double(k);
            y = NormalCDFInverse(y);
            randomNormal(x_vector);
            for (int index = 0; index < 3; index++) {
                z_vector[index] = v[index] * y;
                z_vector[index] += A[index][0] * x_vector[0] + A[index][1] * x_vector[1] + A[index][2] * x_vector[2];
            }
            mean += g_x(z_vector);
        }
    }
    return mean / (num_rep * k);
}

void q2() {
    unsigned long number_of_experiments = 10000;
    unsigned long n = 100;
    unsigned int k = 5;

    std::vector<double> vanillaMCOutput(number_of_experiments);
    std::vector<double> stratifiedMCOutput(number_of_experiments);
    for (int i=0; i < int(number_of_experiments); i++) {
        vanillaMCOutput[i] = vanillaMC2(n);
        stratifiedMCOutput[i] = kFoldStratify(k, n);
    }

    std::cout << "q2:" << std::endl;
    std::cout << "num stratum:\t" << k << std::endl;

    std::cout << "Vanilla MC:\t" << "mean:\t" << mean(vanillaMCOutput) << "\tvariance:\t" << computeVariance(vanillaMCOutput) << std::endl;
    std::cout << "stratified MC:\t" << "mean:\t" << mean(stratifiedMCOutput) << "\tvariance:\t" << computeVariance(stratifiedMCOutput) << std::endl;
    std::cout << "efficiency gain:\t" << computeVariance(vanillaMCOutput) / computeVariance(stratifiedMCOutput) << std::endl;
}