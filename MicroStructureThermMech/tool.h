#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/src/KroneckerProduct/KroneckerTensorProduct.h>
#include <iostream>
#include <vector>
#include <unordered_map>

#include <opencv2/opencv.hpp>

using namespace Eigen;
using namespace std;
using namespace cv;

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L /* pi */
#endif                                               // ! M_PI

MatrixXd homogenize(int lx, int ly, 
    std::vector<double> lambda, 
    std::vector<double> mu, 
    double phi, 
    MatrixXi x);

void elementMatVec(double a,
    double b,
    double phi,
    MatrixXd& keLambda,
    MatrixXd& keMu,
    MatrixXd& feLambda,
    MatrixXd& feMu);

MatrixXi image2matrix(string filename);