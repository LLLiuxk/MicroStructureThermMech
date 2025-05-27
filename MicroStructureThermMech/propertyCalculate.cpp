#include "propertyCalculate.h"
#include <stdexcept>  // std::runtime_error
#include <cmath>      // 示例中使用 sqrt 等
#include <iostream>   // 可用于调试打印日志，例如 std::cout

namespace msGen {

    /**
     * @brief 计算弹性张量示例
     */
    MatrixXd generateCHMatrix(MatrixXi input)
    {
        int lx = input.cols();
        int ly = input.rows();

        double E = 1.0f;
        double mut = 0.45;
        double Lambda = mut * E / ((1 + mut) * (1 - 2 * mut));
        double Mu = E / (2 * (1 + mut));

        std::vector<double> lambda = { Lambda, Lambda * 1e-9 };
        std::vector<double> mu = { Mu, Mu * 1e-9 };
        double phi = 90;
        //MatrixXi x(5, 5);
        //x << 1, 2, 1, 1, 2,
        //	2, 1, 1, 1, 2,
        //	1, 1, 1, 1, 1,
        //	1, 2, 1, 2, 1,
        //	1, 1, 1, 1, 1;
        MatrixXd CH = homogenize(lx, ly, lambda, mu, phi, input);

        std::cout << CH << '\n';
        return CH;
    }

    /**
     * @brief 计算热导率示例
     */
    double calculateThermalConductivity(const std::vector<std::vector<int>>& microstructure)
    {
        // 1. 参数检查
        if (microstructure.empty() || microstructure[0].empty()) {
            throw std::runtime_error("[proCal::calculateThermalConductivity] microstructure 为空。");
        }

        int height = microstructure.size();
        int width = microstructure[0].size();

        // 2. 统计材料“1”所占比
        int count1 = 0;
        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                if (microstructure[row][col] == 1) {
                    count1++;
                }
            }
        }

        double total = static_cast<double>(height * width);
        double ratio1 = static_cast<double>(count1) / total;
        double ratio0 = 1.0 - ratio1;

        // 3. 示例：假设材料“0”热导率 k0 = 0.1, 材料“1”热导率 k1 = 5.0
        //    并做一个简单线性混合。更真实的计算可用并联、串联或更复杂的规则
        double k0 = 0.1;
        double k1 = 5.0;
        // 简单加权，示例用 Voigt 平均 (并联模型)
        // 也可使用 Reuss 平均 (串联模型) 或 Hashin–Shtrikman 界限等
        double k_eff = k0 * ratio0 + k1 * ratio1;

        return k_eff;
    }

} // namespace proCal