#include "Visualization.h"
#include <fstream>     // std::ofstream
#include <algorithm>   // std::min, std::max
#include <stdexcept>   // std::runtime_error
#include <iostream>

namespace Visualization {

    namespace {
        // 辅助函数: 将 val 从区间 [inMin, inMax] 线性映射到 [0, 255]。
        // 若 val 超出 [inMin, inMax]，可做适当裁剪。
        unsigned char linearMapToGray(double val, double inMin, double inMax)
        {
            if (inMax - inMin < 1e-12) {
                // 避免被 0 除错误
                return 0;
            }
            // 裁剪
            if (val < inMin) val = inMin;
            if (val > inMax) val = inMax;

            double ratio = (val - inMin) / (inMax - inMin);
            return static_cast<unsigned char>(ratio * 255.0);
        }

        // 辅助函数: 将一个二维灰度数据保存为 PGM (Portable GrayMap) 格式
        // grayData[row][col] 范围建议在 [0,255]
        // height = grayData.size(), width = grayData[0].size()
        bool saveAsPGM(const std::vector<std::vector<unsigned char>>& grayData,
            const std::string& fileName)
        {
            if (grayData.empty() || grayData[0].empty()) {
                std::cerr << "[saveAsPGM] 输入灰度数据为空，无法保存: " << fileName << std::endl;
                return false;
            }

            int height = static_cast<int>(grayData.size());
            int width = static_cast<int>(grayData[0].size());

            std::ofstream out(fileName, std::ios::binary);
            if (!out.is_open()) {
                std::cerr << "[saveAsPGM] 无法打开文件: " << fileName << std::endl;
                return false;
            }

            // 写 PGM 头部
            // P5 表示 binary PGM
            out << "P5\n" << width << " " << height << "\n255\n";

            // 写像素数据
            for (int row = 0; row < height; ++row) {
                out.write(reinterpret_cast<const char*>(grayData[row].data()), width);
            }

            out.close();
            return true;
        }
    } // end anonymous namespace

    void generateResultImages(
        const std::vector<std::vector<int>>& microstructure,
        const std::vector<double>& elasticTensor,
        double thermalConductivity,
        const std::string& microstructureFileName,
        const std::string& elasticityFileName,
        const std::string& thermalFileName
    )
    {
        //======== 1. 微结构可视化 ========//
        if (microstructure.empty() || microstructure[0].empty()) {
            throw std::runtime_error("[Visualization::generateResultImages] 微结构数据为空。");
        }

        int height = static_cast<int>(microstructure.size());
        int width = static_cast<int>(microstructure[0].size());

        // 将 0/1 转为 [0/255] 灰度，0->0(黑), 1->255(白)
        std::vector<std::vector<unsigned char>> microGray(height, std::vector<unsigned char>(width, 0));
        for (int r = 0; r < height; ++r) {
            for (int c = 0; c < width; ++c) {
                // 这里简单处理: microstructure[r][c] = 0/1
                if (microstructure[r][c] == 1) {
                    microGray[r][c] = 255;  // 白色
                }
                else {
                    microGray[r][c] = 0;    // 黑色
                }
            }
        }

        // 保存为 PGM
        saveAsPGM(microGray, microstructureFileName);
        std::cout << "[Visualization] 已保存微结构图: " << microstructureFileName << std::endl;


        //======== 2. 弹性张量可视化 ========//
        // 假设 elasticTensor 有 4 个分量 (C11, C22, C33, C12)
        // 将它们映射为 2x2 的灰度图
        if (elasticTensor.size() < 4) {
            std::cerr << "[Visualization] 警告: 弹性张量分量数不足 4，无法完整演示 2x2 可视化。" << std::endl;
        }

        int tensorSize = static_cast<int>(elasticTensor.size());
        int nRows = 2;
        int nCols = 2;
        std::vector<std::vector<unsigned char>> tensorGray(nRows, std::vector<unsigned char>(nCols, 0));

        // 找到最大最小值, 用于线性映射
        // 如果分量 < 4，则重复使用前几个或先行判断处理
        double minVal = 1e30;
        double maxVal = -1e30;
        for (int i = 0; i < tensorSize; ++i) {
            if (elasticTensor[i] < minVal) minVal = elasticTensor[i];
            if (elasticTensor[i] > maxVal) maxVal = elasticTensor[i];
        }
        if (tensorSize == 0) {
            // 避免被 0 除错误
            minVal = 0;
            maxVal = 1;
        }

        // 按 (row, col) = (0,0)->C11, (0,1)->C22, (1,0)->C33, (1,1)->C12 示例
        // 如果数量不足 4，就一一对应，对超出的不处理
        std::vector<int> mapIndex = { 0, 1, 2, 3 };
        for (int idx = 0; idx < 4 && idx < tensorSize; ++idx) {
            int r = idx / 2;
            int c = idx % 2;
            double val = elasticTensor[mapIndex[idx]];
            tensorGray[r][c] = linearMapToGray(val, minVal, maxVal);
        }

        // 保存为 PGM
        saveAsPGM(tensorGray, elasticityFileName);
        std::cout << "[Visualization] 已保存弹性张量可视化图: " << elasticityFileName << std::endl;


        //======== 3. 热导率可视化 ========//
        // 这里 thermalConductivity 是单个标量，只需要画一个 64x64 的灰度方块
        // 用 linearMapToGray 做个映射即可
        // 示例：若全部像素都使用同一个灰度值
        std::vector<std::vector<unsigned char>> thermalGray(64, std::vector<unsigned char>(64, 0));

        // 可以先行定义一个你期望的热导率映射范围，示例中简单设 [0, 10] 为映射区间
        double kMin = 0.0;
        double kMax = 10.0;
        unsigned char grayVal = linearMapToGray(thermalConductivity, kMin, kMax);

        for (int r = 0; r < 64; ++r) {
            for (int c = 0; c < 64; ++c) {
                thermalGray[r][c] = grayVal;
            }
        }

        // 保存为 PGM
        saveAsPGM(thermalGray, thermalFileName);
        std::cout << "[Visualization] 已保存热导率可视化图: " << thermalFileName << std::endl;
    }

} // namespace Visualization