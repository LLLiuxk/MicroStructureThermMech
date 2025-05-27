#pragma once
#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <vector>
#include <string>

namespace msGen {

    /**
     * @brief 生成三张图文件，分别展示:
     *  1) 二维黑白微结构
     *  2) 弹性张量可视化图 (示例将弹性张量四个值映射到简单的 2x2 灰度图)
     *  3) 热导率可视化图 (示例用一个尺度图表示单个标量的大小)
     *
     * @param microstructure          二维黑白图 (0,1)
     * @param elasticTensor           一维的弹性张量数据 (示例假定有 4 个分量)
     * @param thermalConductivity     整体热导率 (标量)
     * @param microstructureFileName  微结构可视化保存文件名
     * @param elasticityFileName      弹性张量可视化保存文件名
     * @param thermalFileName         热导率可视化保存文件名
     */
    void generateResultImages(
        const std::vector<std::vector<int>>& microstructure,
        const std::vector<double>& elasticTensor,
        double thermalConductivity,
        const std::string& microstructureFileName,
        const std::string& elasticityFileName,
        const std::string& thermalFileName
    );

} // namespace Visualization

#endif // VISUALIZATION_H