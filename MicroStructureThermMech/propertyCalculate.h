#pragma once
#ifndef PROPERTYCALCULATE_H
#define PROPERTYCALCULATE_H

#include <vector>

#include "tool.h"

using namespace Eigen;

namespace msGen {

    
    /** @brief 根据输入的二维黑白图计算弹性张量
     *
     * @param microstructure 二维黑白图（如 0 表示空或软材料，1 表示实或硬材料）
     * @return 返回一个一维数组，用于表示弹性张量的各个分量。*/
    MatrixXd calculateCH(MatrixXi input, double E = 1.0, double nu = 0.45, double phi = 90);
    MatrixXd calculateCH(string filepath, double E = 1.0, double nu = 0.45, double phi = 90);

    MatrixXd calculateKappaH(MatrixXi input, double mu1 = 10, double mu2 = 1, double phi = 90 );
    MatrixXd calculateKappaH(string filepath, double mu1 = 10, double mu2 = 1, double phi = 90);

} // namespace proCal

#endif // PROPERTYCALCULATE_H