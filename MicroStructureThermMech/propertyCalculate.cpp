#include "propertyCalculate.h"
#include <stdexcept>  // std::runtime_error
#include <cmath>      // ʾ����ʹ�� sqrt ��
#include <iostream>   // �����ڵ��Դ�ӡ��־������ std::cout

namespace msGen {

    /**
     * @brief ���㵯������ʾ��
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
     * @brief �����ȵ���ʾ��
     */
    double calculateThermalConductivity(const std::vector<std::vector<int>>& microstructure)
    {
        // 1. �������
        if (microstructure.empty() || microstructure[0].empty()) {
            throw std::runtime_error("[proCal::calculateThermalConductivity] microstructure Ϊ�ա�");
        }

        int height = microstructure.size();
        int width = microstructure[0].size();

        // 2. ͳ�Ʋ��ϡ�1����ռ��
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

        // 3. ʾ����������ϡ�0���ȵ��� k0 = 0.1, ���ϡ�1���ȵ��� k1 = 5.0
        //    ����һ�������Ի�ϡ�����ʵ�ļ�����ò���������������ӵĹ���
        double k0 = 0.1;
        double k1 = 5.0;
        // �򵥼�Ȩ��ʾ���� Voigt ƽ�� (����ģ��)
        // Ҳ��ʹ�� Reuss ƽ�� (����ģ��) �� Hashin�CShtrikman ���޵�
        double k_eff = k0 * ratio0 + k1 * ratio1;

        return k_eff;
    }

} // namespace proCal