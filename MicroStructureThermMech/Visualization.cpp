#include "Visualization.h"
#include <fstream>     // std::ofstream
#include <algorithm>   // std::min, std::max
#include <stdexcept>   // std::runtime_error
#include <iostream>

namespace Visualization {

    namespace {
        // ��������: �� val ������ [inMin, inMax] ����ӳ�䵽 [0, 255]��
        // �� val ���� [inMin, inMax]�������ʵ��ü���
        unsigned char linearMapToGray(double val, double inMin, double inMax)
        {
            if (inMax - inMin < 1e-12) {
                // ���ⱻ 0 ������
                return 0;
            }
            // �ü�
            if (val < inMin) val = inMin;
            if (val > inMax) val = inMax;

            double ratio = (val - inMin) / (inMax - inMin);
            return static_cast<unsigned char>(ratio * 255.0);
        }

        // ��������: ��һ����ά�Ҷ����ݱ���Ϊ PGM (Portable GrayMap) ��ʽ
        // grayData[row][col] ��Χ������ [0,255]
        // height = grayData.size(), width = grayData[0].size()
        bool saveAsPGM(const std::vector<std::vector<unsigned char>>& grayData,
            const std::string& fileName)
        {
            if (grayData.empty() || grayData[0].empty()) {
                std::cerr << "[saveAsPGM] ����Ҷ�����Ϊ�գ��޷�����: " << fileName << std::endl;
                return false;
            }

            int height = static_cast<int>(grayData.size());
            int width = static_cast<int>(grayData[0].size());

            std::ofstream out(fileName, std::ios::binary);
            if (!out.is_open()) {
                std::cerr << "[saveAsPGM] �޷����ļ�: " << fileName << std::endl;
                return false;
            }

            // д PGM ͷ��
            // P5 ��ʾ binary PGM
            out << "P5\n" << width << " " << height << "\n255\n";

            // д��������
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
        //======== 1. ΢�ṹ���ӻ� ========//
        if (microstructure.empty() || microstructure[0].empty()) {
            throw std::runtime_error("[Visualization::generateResultImages] ΢�ṹ����Ϊ�ա�");
        }

        int height = static_cast<int>(microstructure.size());
        int width = static_cast<int>(microstructure[0].size());

        // �� 0/1 תΪ [0/255] �Ҷȣ�0->0(��), 1->255(��)
        std::vector<std::vector<unsigned char>> microGray(height, std::vector<unsigned char>(width, 0));
        for (int r = 0; r < height; ++r) {
            for (int c = 0; c < width; ++c) {
                // ����򵥴���: microstructure[r][c] = 0/1
                if (microstructure[r][c] == 1) {
                    microGray[r][c] = 255;  // ��ɫ
                }
                else {
                    microGray[r][c] = 0;    // ��ɫ
                }
            }
        }

        // ����Ϊ PGM
        saveAsPGM(microGray, microstructureFileName);
        std::cout << "[Visualization] �ѱ���΢�ṹͼ: " << microstructureFileName << std::endl;


        //======== 2. �����������ӻ� ========//
        // ���� elasticTensor �� 4 ������ (C11, C22, C33, C12)
        // ������ӳ��Ϊ 2x2 �ĻҶ�ͼ
        if (elasticTensor.size() < 4) {
            std::cerr << "[Visualization] ����: ������������������ 4���޷�������ʾ 2x2 ���ӻ���" << std::endl;
        }

        int tensorSize = static_cast<int>(elasticTensor.size());
        int nRows = 2;
        int nCols = 2;
        std::vector<std::vector<unsigned char>> tensorGray(nRows, std::vector<unsigned char>(nCols, 0));

        // �ҵ������Сֵ, ��������ӳ��
        // ������� < 4�����ظ�ʹ��ǰ�����������жϴ���
        double minVal = 1e30;
        double maxVal = -1e30;
        for (int i = 0; i < tensorSize; ++i) {
            if (elasticTensor[i] < minVal) minVal = elasticTensor[i];
            if (elasticTensor[i] > maxVal) maxVal = elasticTensor[i];
        }
        if (tensorSize == 0) {
            // ���ⱻ 0 ������
            minVal = 0;
            maxVal = 1;
        }

        // �� (row, col) = (0,0)->C11, (0,1)->C22, (1,0)->C33, (1,1)->C12 ʾ��
        // ����������� 4����һһ��Ӧ���Գ����Ĳ�����
        std::vector<int> mapIndex = { 0, 1, 2, 3 };
        for (int idx = 0; idx < 4 && idx < tensorSize; ++idx) {
            int r = idx / 2;
            int c = idx % 2;
            double val = elasticTensor[mapIndex[idx]];
            tensorGray[r][c] = linearMapToGray(val, minVal, maxVal);
        }

        // ����Ϊ PGM
        saveAsPGM(tensorGray, elasticityFileName);
        std::cout << "[Visualization] �ѱ��浯���������ӻ�ͼ: " << elasticityFileName << std::endl;


        //======== 3. �ȵ��ʿ��ӻ� ========//
        // ���� thermalConductivity �ǵ���������ֻ��Ҫ��һ�� 64x64 �ĻҶȷ���
        // �� linearMapToGray ����ӳ�伴��
        // ʾ������ȫ�����ض�ʹ��ͬһ���Ҷ�ֵ
        std::vector<std::vector<unsigned char>> thermalGray(64, std::vector<unsigned char>(64, 0));

        // �������ж���һ�����������ȵ���ӳ�䷶Χ��ʾ���м��� [0, 10] Ϊӳ������
        double kMin = 0.0;
        double kMax = 10.0;
        unsigned char grayVal = linearMapToGray(thermalConductivity, kMin, kMax);

        for (int r = 0; r < 64; ++r) {
            for (int c = 0; c < 64; ++c) {
                thermalGray[r][c] = grayVal;
            }
        }

        // ����Ϊ PGM
        saveAsPGM(thermalGray, thermalFileName);
        std::cout << "[Visualization] �ѱ����ȵ��ʿ��ӻ�ͼ: " << thermalFileName << std::endl;
    }

} // namespace Visualization