#include "propertyCalculate.h"
#include <stdexcept>  // std::runtime_error
#include <cmath>      // 示例中使用 sqrt 等
#include <iostream>   // 可用于调试打印日志，例如 std::cout

namespace msGen {

    MatrixXd calculateCH(MatrixXi input, double E, double nu, double phi )
    {
        int lx = input.cols();
        int ly = input.rows();

        double lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
        double mu = E / (2 * (1 + nu));

        std::vector<double> Lambda = { lambda, lambda * 1e-9 };
        std::vector<double> Mu = { mu, mu * 1e-9 };
        
        MatrixXd CH = homogenize(lx, ly, Lambda, Mu, phi, input);

        //std::cout << "CH:"<<CH << std::endl;
        return CH;
    }

    MatrixXd calculateCH(string filepath, double E, double nu, double phi)
    {
        MatrixXi binaryMat = image2matrix(filepath);

        std::vector<double> Lambda = { nu * E / ((1 + nu) * (1 - 2 * nu)) , 1e-20 };
        std::vector<double> Mu = { E / (2 * (1 + nu)) , 1e-20 };

        MatrixXd CH = homogenize(binaryMat.cols(), binaryMat.rows(), Lambda, Mu, phi, binaryMat);
        
        //std::cout << "CH:" << CH << std::endl;
        return CH;
    }

    MatrixXd calculateKappaH(MatrixXi input, double mu1 , double mu2, double phi )
    {
        int lx = input.cols();
        int ly = input.rows();
        std::vector<double> Lambda2 = { 0, 0 };
        std::vector<double> Mu2 = { mu1, mu2 };
        MatrixXd KH = homogenize_therm(lx, ly, Lambda2, Mu2, phi, input);
        //cout << "KH:" << endl << KH << endl; 
        return KH;
    }

    MatrixXd calculateKappaH(string filepath, double mu1, double mu2 , double phi )
    {
        Eigen::MatrixXi binaryMat = image2matrix(filepath);
        std::vector<double> Lambda2 = { 0, 0 };
        std::vector<double> Mu2 = { mu1, mu2 };
        MatrixXd KH = homogenize_therm(binaryMat.cols(), binaryMat.rows(), Lambda2, Mu2, phi, binaryMat);
       // cout << "KH:" << endl << KH << endl;
        return KH;
    }
    
} // namespace msGen