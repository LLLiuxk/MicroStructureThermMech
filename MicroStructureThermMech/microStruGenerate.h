#pragma once
#ifndef MICROSTRUGENERATE_H
#define MICROSTRUGENERATE_H

#include <vector>

enum ConnectionMode {
    MODE_NONE = 0,
    MODE_LINEAR,
    MODE_BSPLINE,
    MODE_RANDOM,
    // ... 其他模式
};

namespace msGen {

    /**
     * @brief 生成微结构的函数示例。
     *
     * 该函数示例使用 8 个参数来表示微结构的关键要素：
     * 1) leftPointCount         左边分布点的数量
     * 2) topPointCount          上边分布点的数量
     * 3) leftPositions          左边每个点在垂直方向上的位置（假设 y 方向）
     * 4) topPositions           上边每个点在水平方向上的位置（假设 x 方向）
     * 5) connectionWidthLeft    从左边出发的连接宽度
     * 6) connectionWidthTop     从上边出发的连接宽度
     * 7) connectionModeLeft     左边分布点之间的连接方式（可以自定义枚举或整型）
     * 8) connectionModeTop      上边分布点之间的连接方式（可以自定义枚举或整型）
     *
     * @param leftPointCount         左边分布点数量
     * @param topPointCount          上边分布点数量
     * @param leftPositions          左边分布点的位置向量
     * @param topPositions           上边分布点的位置向量
     * @param connectionWidthLeft    左侧连接线宽度（可理解为画线时的笔触宽度等）
     * @param connectionWidthTop     上侧连接线宽度
     * @param connectionModeLeft     左侧点之间的连接方式（可以用枚举或 int 表示不同的连接模式）
     * @param connectionModeTop      上侧点之间的连接方式
     *
     * @return 以二维整型数组(如 0/1)返回生成后的微结构示意图,
     *         0 可能表示空白区域，1 可能表示材料区域。
     *         可根据实际项目需求改为 double 或其他容器类型。
     */
    std::vector<std::vector<int>> generateMicrostructure(
        int leftPointCount,
        int topPointCount,
        const std::vector<double>& leftPositions,
        const std::vector<double>& topPositions,
        double connectionWidthLeft,
        double connectionWidthTop,
        int connectionModeLeft,
        int connectionModeTop
    );

} // namespace msGen

#endif // MICROSTRUGENERATE_H