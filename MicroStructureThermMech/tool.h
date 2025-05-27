#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/src/KroneckerProduct/KroneckerTensorProduct.h>
#include <iostream>
#include <vector>
using namespace Eigen;

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L /* pi */
#endif                                               // ! M_PI

MatrixXd homogenize(int lx, int ly, 
    std::vector<double> lambda, 
    std::vector<double> mu, 
    double phi, 
    MatrixXi x);

MatrixXd generateCHMatrix(MatrixXi input);

void elementMatVec(double a,
    double b,
    double phi,
    MatrixXd& keLambda,
    MatrixXd& keMu,
    MatrixXd& feLambda,
    MatrixXd& feMu);