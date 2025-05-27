#pragma once
#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <vector>
#include <string>

namespace msGen {

    /**
     * @brief ��������ͼ�ļ����ֱ�չʾ:
     *  1) ��ά�ڰ�΢�ṹ
     *  2) �����������ӻ�ͼ (ʾ�������������ĸ�ֵӳ�䵽�򵥵� 2x2 �Ҷ�ͼ)
     *  3) �ȵ��ʿ��ӻ�ͼ (ʾ����һ���߶�ͼ��ʾ���������Ĵ�С)
     *
     * @param microstructure          ��ά�ڰ�ͼ (0,1)
     * @param elasticTensor           һά�ĵ����������� (ʾ���ٶ��� 4 ������)
     * @param thermalConductivity     �����ȵ��� (����)
     * @param microstructureFileName  ΢�ṹ���ӻ������ļ���
     * @param elasticityFileName      �����������ӻ������ļ���
     * @param thermalFileName         �ȵ��ʿ��ӻ������ļ���
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