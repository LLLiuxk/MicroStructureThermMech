#include <iostream>
#include <vector>
#include <string>


#include "tool.h"
#include "microStruGenerate.h"
#include "propertyCalculate.h"
#include "Visualization.h"

int main()
{
    // 设置微结构设计所需的边界参数
    int width = 100;              // 微结构图的宽度
    int height = 100;             // 微结构图的高度
    double porosity = 0.3;        //孔隙率
    int seed = 42;                // 随机种子，用于生成微结构时的可重复性
    int leftPointCount;
    int topPointCount;
    const std::vector<double> leftPositions;
    const std::vector<double> topPositions;
    double connectionWidthLeft;
    double connectionWidthTop;
    int connectionModeLeft;
    int connectionModeTop;

    // 生成黑白二维图,返回一个二维容器，如std::vector<std::vector<int>>， 其中 0 代表黑色，1 代表白色
    auto microstructure = msGen::generateMicrostructure(
        leftPointCount, topPointCount, leftPositions, topPositions, connectionWidthLeft, connectionWidthTop, connectionModeLeft,  connectionModeTop
    );

    // 计算每个二维微结构的热导率和弹性张量
    double thermalConductivity = proCal::computeThermalConductivity(
        microstructure
    );
    auto elasticTensor = proCal::computeElasticTensor(
        microstructure
    );

    // 可视化,这里分别可视化 (1) 黑白微结构图；(2) 热导率分布图；(3) 弹性张量分布图
   
    Visualization::plotMicrostructure(
        microstructure, "Microstructure.png"
    );


    std::cout << "热导率计算结果: " << thermalConductivity << std::endl;
    std::cout << "弹性张量主元: ";
    for (auto val : elasticTensor) {
        std::cout << val << " ";
    }

    std::cout << std::endl<<"微结构设计、性质计算与可视化全部完成！" << std::endl;

    return 0;
}

