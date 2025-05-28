#include "microStruGenerate.h"

#include <algorithm>  // std::max_element
#include <cmath>      // std::round, std::abs
#include <stdexcept>  // std::runtime_error
#include <iostream>   // ����ʱ����Դ�ӡ��־

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
         * @brief ��ȡ����������ʾ������һά�����洢�����Զ���Ϊ�����������ʽ����
         * @return ���� m_elasticTensor������ [C11, C22, C33, C12]��
         */
        const std::vector<double>& getElasticTensor() const;

        /**
         * @brief ��ȡ�ȵ���������ʾ������һά�����洢������ 2��2����
         * @return ���� m_thermalConductivityTensor������ [k11, k22, k12, k21]��
         */
        const std::vector<double>& getThermalConductivityTensor() const;

    private:
        // ΢�ṹ����
        int m_leftPointCount;
        int m_topPointCount;
        std::vector<double> m_leftPositions;
        std::vector<double> m_topPositions;
        double m_connectionWidthLeft;
        double m_connectionWidthTop;
        int m_connectionModeLeft;
        int m_connectionModeTop;

        // ���ɵĶ�ά�ڰ�ͼ (0/1)
        std::vector<std::vector<int>> m_microstructure;

        // ����������ʾ����4 ������ [C11, C22, C33, C12]��
        std::vector<double> m_elasticTensor;

        // �ȵ���������ʾ����2��2���� 4 ������ [k11, k22, k12, k21]��
        std::vector<double> m_thermalConductivityTensor;

    private:
        /**
         * @brief ���㵯������ (ʾ�������滻Ϊʵ�ʵ���ѧ���ʻ��㷨)
         */
        void computeElasticTensor();

        /**
         * @brief �����ȵ������� (ʾ�������滻Ϊʵ�ʵ��ȴ������ʻ��㷨)
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