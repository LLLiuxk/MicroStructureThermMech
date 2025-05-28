#include "microStruGenerate.h"

#include <algorithm>  // std::max_element
#include <cmath>      // std::round, std::abs
#include <stdexcept>  // std::runtime_error
#include <iostream>   // 调试时候可以打印日志

namespace msGen {

    class Microstructure
    {
    public:
      
        Microstructure(
            int leftPointCount,
            int topPointCount,
            const std::vector<double>& leftPositions,
            const std::vector<double>& topPositions,
            double connectionWidthLeft,
            double connectionWidthTop,
            int connectionModeLeft,
            int connectionModeTop
        );

        ~Microstructure() {};
        
        void generate();

        
        const std::vector<std::vector<int>>& getMicrostructure() const;

        /**
         * @brief 获取弹性张量（示例中以一维向量存储，可自定义为矩阵或其他格式）。
         * @return 返回 m_elasticTensor，例如 [C11, C22, C33, C12]。
         */
        const std::vector<double>& getElasticTensor() const;

        /**
         * @brief 获取热导率张量（示例中以一维向量存储，假设 2×2）。
         * @return 返回 m_thermalConductivityTensor，例如 [k11, k22, k12, k21]。
         */
        const std::vector<double>& getThermalConductivityTensor() const;

    private:
        // 微结构参数
        int m_leftPointCount;
        int m_topPointCount;
        std::vector<double> m_leftPositions;
        std::vector<double> m_topPositions;
        double m_connectionWidthLeft;
        double m_connectionWidthTop;
        int m_connectionModeLeft;
        int m_connectionModeTop;

        // 生成的二维黑白图 (0/1)
        std::vector<std::vector<int>> m_microstructure;

        // 弹性张量（示例：4 个分量 [C11, C22, C33, C12]）
        std::vector<double> m_elasticTensor;

        // 热导率张量（示例：2×2，共 4 个分量 [k11, k22, k12, k21]）
        std::vector<double> m_thermalConductivityTensor;

    private:
        /**
         * @brief 计算弹性张量 (示例，可替换为实质的力学均质化算法)
         */
        void computeElasticTensor();

        /**
         * @brief 计算热导率张量 (示例，可替换为实质的热传导均质化算法)
         */
        void computeThermalConductivityTensor();

    };




    std::vector<std::vector<int>> generateMicrostructure(
        int leftPointCount,
        int topPointCount,
        const std::vector<double>& leftPositions,
        const std::vector<double>& topPositions,
        double connectionWidthLeft,
        double connectionWidthTop,
        int connectionModeLeft,
        int connectionModeTop
    )
    {
        std::vector<std::vector<int>> microstructure;
        return microstructure;
    }



} 