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
    MatrixXd generateCHMatrix(MatrixXi input);

    /**
     * @brief 根据输入的二维黑白图计算热导率
     *
     * @param microstructure 二维黑白图
     * @return 返回一个双精度浮点数，表示整体热导率值
     */
    double calculateThermalConductivity(const std::vector<std::vector<int>>& microstructure);

} // namespace proCal

#endif // PROPERTYCALCULATE_H