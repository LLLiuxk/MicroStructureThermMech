#pragma once
#ifndef MICROSTRUGENERATE_H
#define MICROSTRUGENERATE_H

#include <vector>

enum ConnectionMode {
    MODE_NONE = 0,
    MODE_LINEAR,
    MODE_BSPLINE,
    MODE_RANDOM,
    // ... ����ģʽ
};

namespace msGen {

    /**
     * @brief ����΢�ṹ�ĺ���ʾ����
     *
     * �ú���ʾ��ʹ�� 8 ����������ʾ΢�ṹ�Ĺؼ�Ҫ�أ�
     * 1) leftPointCount         ��߷ֲ��������
     * 2) topPointCount          �ϱ߷ֲ��������
     * 3) leftPositions          ���ÿ�����ڴ�ֱ�����ϵ�λ�ã����� y ����
     * 4) topPositions           �ϱ�ÿ������ˮƽ�����ϵ�λ�ã����� x ����
     * 5) connectionWidthLeft    ����߳��������ӿ��
     * 6) connectionWidthTop     ���ϱ߳��������ӿ��
     * 7) connectionModeLeft     ��߷ֲ���֮������ӷ�ʽ�������Զ���ö�ٻ����ͣ�
     * 8) connectionModeTop      �ϱ߷ֲ���֮������ӷ�ʽ�������Զ���ö�ٻ����ͣ�
     *
     * @param leftPointCount         ��߷ֲ�������
     * @param topPointCount          �ϱ߷ֲ�������
     * @param leftPositions          ��߷ֲ����λ������
     * @param topPositions           �ϱ߷ֲ����λ������
     * @param connectionWidthLeft    ��������߿�ȣ������Ϊ����ʱ�ıʴ���ȵȣ�
     * @param connectionWidthTop     �ϲ������߿��
     * @param connectionModeLeft     ����֮������ӷ�ʽ��������ö�ٻ� int ��ʾ��ͬ������ģʽ��
     * @param connectionModeTop      �ϲ��֮������ӷ�ʽ
     *
     * @return �Զ�ά��������(�� 0/1)�������ɺ��΢�ṹʾ��ͼ,
     *         0 ���ܱ�ʾ�հ�����1 ���ܱ�ʾ��������
     *         �ɸ���ʵ����Ŀ�����Ϊ double �������������͡�
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