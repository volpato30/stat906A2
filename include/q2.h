//
// Created by rui on 28/10/18.
//

#ifndef STAT906A2_Q2_H
#define STAT906A2_Q2_H

#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <random>
// copy the implementation of CDF
// Inverse from https://www.johndcook.com/blog/cpp_phi_inverse/

const double v[3] = {0.63245553,  0.63245553, -0.4472136};
const double A[3][3] = {
        {0.6, -0.4, 0.28284271},
        {-0.4, 0.6, 0.28284271},
        {0.2828471, 0.28284271, 0.8}
};


void q2();

#endif //STAT906A2_Q2_H
