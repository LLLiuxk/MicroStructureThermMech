#pragma once
#ifndef PROPERTYCALCULATE_H
#define PROPERTYCALCULATE_H

#include <vector>

#include "tool.h"

using namespace Eigen;

namespace msGen {

    
    /** @brief ��������Ķ�ά�ڰ�ͼ���㵯������
     *
     * @param microstructure ��ά�ڰ�ͼ���� 0 ��ʾ�ջ�����ϣ�1 ��ʾʵ��Ӳ���ϣ�
     * @return ����һ��һά���飬���ڱ�ʾ���������ĸ���������*/
    MatrixXd generateCHMatrix(MatrixXi input);

    /**
     * @brief ��������Ķ�ά�ڰ�ͼ�����ȵ���
     *
     * @param microstructure ��ά�ڰ�ͼ
     * @return ����һ��˫���ȸ���������ʾ�����ȵ���ֵ
     */
    double calculateThermalConductivity(const std::vector<std::vector<int>>& microstructure);

} // namespace proCal

#endif // PROPERTYCALCULATE_H