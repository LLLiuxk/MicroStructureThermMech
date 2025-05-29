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
    MatrixXd calculateCH(MatrixXi input, double E = 1.0, double nu = 0.45, double phi = 90);
    MatrixXd calculateCH(string filepath, double E = 1.0, double nu = 0.45, double phi = 90);

    MatrixXd calculateKappaH(MatrixXi input, double mu1 = 10, double mu2 = 1, double phi = 90 );
    MatrixXd calculateKappaH(string filepath, double mu1 = 10, double mu2 = 1, double phi = 90);

} // namespace proCal

#endif // PROPERTYCALCULATE_H