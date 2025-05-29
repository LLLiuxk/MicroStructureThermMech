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

// calculate homo property
MatrixXd homogenize(int lx, int ly, 
    std::vector<double> lambda, 
    std::vector<double> mu, 
    double phi, 
    MatrixXi x);

MatrixXd homogenize_therm(int lx, int ly,
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

void elementMatVec_therm(double a,
    double b,
    double phi,
    MatrixXd& keLambda,
    MatrixXd& keMu,
    MatrixXd& feLambda,
    MatrixXd& feMu);

MatrixXi image2matrix(string filename);

void showSparseMatrix(SparseMatrix<double> X);

//math tools
class BSampleFunction
{
public:
    BSampleFunction() {};
    std::vector<Eigen::Vector2d> ThreeOrderBSplineInterpolatePt(std::vector<Eigen::Vector2d>& pt, int InsertNum);

private:
    double F03(double t);
    double F13(double t);
    double F23(double t);
    double F33(double t);
};

class HermitFunction
{
public:
    HermitFunction() {};

    HermitFunction(const Eigen::MatrixXd& control_points, const Eigen::MatrixXd& tangents);

    HermitFunction(Eigen::Vector2d b, Eigen::Vector2d p0, Eigen::Vector2d e, Eigen::Vector2d p1);

    HermitFunction(Eigen::Vector2d b, Eigen::Vector2d e, double angleb, double anglee);

    void HermitFunction2(Eigen::Vector2d b, Eigen::Vector2d e, double angleb, double anglee);

    void UsingpointWithtangent(Eigen::Vector2d b, Eigen::Vector2d p0, Eigen::Vector2d e, Eigen::Vector2d p1);

    Eigen::Vector2d getpoint(double t);

    Eigen::Vector2d getderivation(double t);
    Eigen::Vector2d GetPoint(double t);

    std::vector<Eigen::Vector2d> getPointsvec(int Numpoints);

    std::vector<Eigen::Vector2d> getPointsvec(std::vector<Eigen::Vector2d> points, std::vector<Eigen::Vector2d> tangents, int Numpoints);

    double distance(const Eigen::Vector2d& point);

private:
    Eigen::MatrixXd m_P; // control points
    Eigen::MatrixXd m_T; // tangents
};