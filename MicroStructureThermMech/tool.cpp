#include "tool.h"

MatrixXi image2matrix(string filename)
{
    cv::Mat grayImage = cv::imread(filename, cv::IMREAD_GRAYSCALE);
    cv::imshow("grayImage", grayImage);
    
    if (grayImage.empty()) {
        std::cerr << "Error: Could not read image " << filename << std::endl;
        return Eigen::MatrixXi(0, 0);
    }

    // 二值化处理 (阈值150，对应MATLAB的im2bw逻辑)
    cv::Mat binaryImage;
    cv::threshold(grayImage, binaryImage, 149, 1, cv::THRESH_BINARY); // <150 -> 0, >=150 -> 1

    // 按MATLAB代码要求 +1（可选）
    binaryImage += 1;  // 将0/1转为1/2

    // 转换为Eigen MatrixXi
    Eigen::MatrixXi result(binaryImage.rows, binaryImage.cols);
    for (int i = 0; i < binaryImage.rows; ++i) {
        for (int j = 0; j < binaryImage.cols; ++j) {
            result(i, j) = binaryImage.at<uchar>(i, j);
        }
    }
    return result;
}


MatrixXd homogenize(
    int lx, int ly, std::vector<double> lambda, std::vector<double> mu, double phi, MatrixXi x)
{
    //INITIALIZE
    // Deduce discretization
    int nely = x.rows(), nelx = x.cols();
    // Stiffness matrix consists of two parts, one belonging to lambda and
    // one belonging to mu. Same goes for load vector
    double dx = (double)lx / nelx, dy = (double)ly / nely;
    int nel = nelx * nely;
    MatrixXd keLambda, keMu, feLambda, feMu;
    elementMatVec(dx / 2, dy / 2, phi, keLambda, keMu, feLambda, feMu);

    // Node numbers and element degrees of freedom for full (not periodic) mesh
    //nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
    MatrixXi nodenrs(nely + 1, nelx + 1);
    for (int i = 0, c = 1; i <= nelx; ++i) {
        for (int j = 0; j <= nely; ++j) {
            nodenrs(j, i) = c++;
        }
    }
    //edofVec = reshape(2 * nodenrs(1:end - 1, 1 : end - 1) + 1, nel, 1);
    VectorXi edofVec(nel);
    for (int i = 0, c = 0; i < nelx; ++i) {
        for (int j = 0; j < nely; ++j) {
            edofVec(c) = 2 * nodenrs(j, i) + 1;
            c++;
        }
    }
    //edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2 * nely + [2 3 0 1] - 2 - 1], nel, 1);
    MatrixXi edofMat(nel, 8);
    edofMat.col(0) = edofVec;
    edofMat.col(1) = edofVec.array() + 1;
    edofMat.col(2) = edofVec.array() + 2 * nely + 2;
    edofMat.col(3) = edofVec.array() + 2 * nely + 3;
    edofMat.col(4) = edofVec.array() + 2 * nely;
    edofMat.col(5) = edofVec.array() + 2 * nely + 1;
    edofMat.col(6) = edofVec.array() - 2;
    edofMat.col(7) = edofVec.array() - 1;
    //cout << edofMat << endl;
    //--------------------------Test----------------------------------------
    //std::cout << "keLambda:\n" << keLambda << "++++++++++++++++++++++\n";
    //std::cout << "keMu:\n" << keMu << "++++++++++++++++++++++\n";
    //std::cout << "feLambda:\n" << feLambda << "++++++++++++++++++++++\n";
    //std::cout << "feMu:\n" << feMu << "++++++++++++++++++++++\n";

    //std::cout << "nodenrs:\n" << nodenrs << "++++++++++++++++++++++\n";
    //std::cout << "edofVec:\n" << edofVec << "++++++++++++++++++++++\n";
    //std::cout << "edofMat:\n" << edofMat << "++++++++++++++++++++++\n";
    //----------------------------------------------------------------------
    //IMPOSE PERIODIC BOUNDARY CONDITIONS
    //Use original edofMat to index into list with the periodic dofs
    int nn = (nelx + 1) * (nely + 1); // Total number of nodes
    int nnP = nelx * nely;            // Total number of unique nodes
    Eigen::MatrixXi nnPArray(nely + 1, nelx + 1);
    nnPArray.setZero();
    for (int i = 0; i < nely; i++)
        for (int j = 0; j < nelx; j++)
            nnPArray(j, i) = i * nelx + j + 1;
    // Extend with a mirror of the top border
    for (int j = 0; j < nelx + 1; j++) {
        nnPArray(nely, j) = nnPArray(0, j);
    }
    // Extend with a mirror of the left border
    for (int i = 0; i < nely + 1; i++) {
        nnPArray(i, nelx) = nnPArray(i, 0);
    }

    //cout << "nnPArray:" << nnPArray << endl;
    Eigen::VectorXd nnParrayVector(nn);
    int indexOfnnParrayVector = 0;
    for (int i = 0; i <= nelx; i++) {
        for (int j = 0; j <= nely; j++) {
            nnParrayVector(indexOfnnParrayVector++) = nnPArray(j, i);
        }
    }

    Eigen::VectorXd dofVector(2 * nn);
    dofVector.setZero();

    for (int i = 0; i < nn; i++) {
        dofVector(i * 2) = 2 * nnParrayVector(i) - 1;
        dofVector(i * 2 + 1) = 2 * nnParrayVector(i);
    }

    Eigen::MatrixXi edofMatNew(nel, 8);
    for (int i = 0; i < nel; i++)
        for (int j = 0; j < 8; j++)
            edofMatNew(i, j) = dofVector(edofMat(i, j) - 1);
    edofMat = edofMatNew;
    int ndof = 2 * nnP; // Number of dofs
//    // ASSEMBLE STIFFNESS MATRIX
//    // Indexing vectors

    Eigen::MatrixXi iKones = Eigen::MatrixXi::Ones(8, 1);
    Eigen::MatrixXi iK = Eigen::kroneckerProduct(edofMat, iKones).transpose();

    Eigen::MatrixXi jKones = Eigen::MatrixXi::Ones(1, 8);
    Eigen::MatrixXi jK = Eigen::kroneckerProduct(edofMat, jKones).transpose();

    //std::cout << "iK:\n" << iK << "++++++++++++++++++++++\n";
    //std::cout << "jK:\n" << jK << "++++++++++++++++++++++\n";

    // Material properties in the different elements
    Eigen::VectorXd lambdaOfx(nel);
    Eigen::VectorXd muOfx(nel);
    for (int i = 0; i < nel; i++) {
        if (x(i) == 1) {
            lambdaOfx(i) = lambda[0];
            muOfx(i) = mu[0];
        }
        else if (x(i) == 2) {
            lambdaOfx(i) = lambda[1];
            muOfx(i) = mu[1];
        }
    }
    // The corresponding stiffness matrix entries
    Eigen::VectorXd keLambdaVector(keLambda.rows() * keLambda.cols());
    Eigen::VectorXd keMuVector(keMu.rows() * keMu.cols());
    int indexOfkeLambdaVector = 0;
    for (int i = 0; i < keLambda.cols(); i++) {
        for (int j = 0; j < keLambda.rows(); j++) {
            keLambdaVector(indexOfkeLambdaVector++) = keLambda(j, i);
        }
    }
    int indexOfkeMuVector = 0;
    for (int i = 0; i < keMu.cols(); i++) {
        for (int j = 0; j < keMu.rows(); j++) {
            keMuVector(indexOfkeMuVector++) = keMu(j, i);
        }
    }
    Eigen::MatrixXd sK = (keLambdaVector * lambdaOfx.transpose() + keMuVector * muOfx.transpose());

    //std::cout << "sK:\n" << sK << "++++++++++++++++++++++\n";

    //std::cout << iK.rows() * iK.cols() << "   ===    " << jK.rows() * jK.cols() << "  ===  " << sK.rows() * sK.cols() << '\n';
    //Map<RowVectorXf> v1(M1.data(), M1.size());
    //cout << "v1:" << endl << v1 << endl;

    //Matrix<float, Dynamic, Dynamic, RowMajor> M2(M1);
    //Map<RowVectorXf> v2(M2.data(), M2.size());
    //cout << "v2:" << endl << v2 << endl;

    Map<VectorXi> iKVector(iK.data(), iK.size());
    Map<VectorXi> jKVector(jK.data(), jK.size());
    Map<VectorXd> sKVector(sK.data(), sK.size());

    Eigen::SparseMatrix<double> K = Eigen::SparseMatrix<double>(ndof, ndof);
    //// Fill K with triplets
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < iK.size(); ++i) {
        triplets.emplace_back(iKVector(i) - 1, jKVector(i) - 1, sKVector(i));
    }
    K.setFromTriplets(triplets.begin(), triplets.end());

    //--------------------------Test----------------------------------------
    //----------------------------------------------------------------------

    //std::cout << "K:\n" << MatrixXd(K) << "++++++++++++++++++++++\n";

    // 需要转置吗？

    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //LOAD VECTORS AND SOLUTION
    Eigen::VectorXd feLambdaVector(feLambda.rows() * feLambda.cols());
    Eigen::VectorXd feMuVector(feMu.rows() * feMu.cols());
    int indexOffeLambdaVector = 0;
    for (int i = 0; i < feLambda.cols(); i++) {
        for (int j = 0; j < feLambda.rows(); j++) {
            feLambdaVector(indexOffeLambdaVector++) = feLambda(j, i);
        }
    }
    int indexOffeMuVector = 0;
    for (int i = 0; i < feMu.cols(); i++) {
        for (int j = 0; j < feMu.rows(); j++) {
            feMuVector(indexOffeMuVector++) = feMu(j, i);
        }
    }

    Eigen::MatrixXd sF = (feLambdaVector * lambdaOfx.transpose() + feMuVector * muOfx.transpose());
    //std::cout << "feLambdaVector:\n" << feLambdaVector << "++++++++++++++++++++++\n";
    //std::cout << "feMuVector:\n" << feMuVector << "++++++++++++++++++++++\n";
    //std::cout << "sF:\n" << sF << "++++++++++++++++++++++\n";

    			//nel * 8
    Eigen::MatrixXi iF(24, nel);
    Eigen::MatrixXi edofMatT = edofMat.transpose();
    iF.block(0, 0, edofMatT.rows(), edofMatT.cols()) = edofMatT;
    iF.block(edofMatT.rows(), 0, edofMatT.rows(), edofMatT.cols()) = edofMatT;
    iF.block(edofMatT.rows() * 2, 0, edofMatT.rows(), edofMatT.cols()) = edofMatT;

    Eigen::MatrixXi jF(24, nel);
    Eigen::MatrixXi jFzeros = Eigen::MatrixXi::Ones(8, nel);
    jF.block(0, 0, jFzeros.rows(), jFzeros.cols()) = jFzeros;
    jF.block(jFzeros.rows(), 0, jFzeros.rows(), jFzeros.cols()) = jFzeros * 2;
    jF.block(jFzeros.rows() * 2, 0, jFzeros.rows(), jFzeros.cols()) = jFzeros * 3;

    Map<VectorXi> iFVector(iF.data(), iF.size());
    Map<VectorXi> jFVector(jF.data(), jF.size());
    Map<VectorXd> sFVector(sF.data(), sF.size());

    Eigen::SparseMatrix<double> F(ndof, 3);

    std::vector<Eigen::Triplet<double>> triplets2;
    for (int i = 0; i < iF.size(); ++i) {
        triplets2.emplace_back(iFVector(i) - 1, jFVector(i) - 1, sFVector(i));
    }
    F.setFromTriplets(triplets2.begin(), triplets2.end());

    SparseLU<SparseMatrix<double>> solver;
    solver.compute(K.block(2, 2, ndof - 2, ndof - 2)); // 求解K中3到ndof的部分，注意要去除边界上的自由度
    SparseMatrix<double> chi(ndof, 3);
    chi = solver.solve(F.block(2, 0, ndof - 2, 3)); // 求解线性方程组

    //showSparseMatrix(K.transpose().block(2, 2, ndof - 2, ndof - 2));
    //showSparseMatrix(X2);
//HOMOGENIZATION
//注意问题  matlab的下标可能与Eigen不同

    Eigen::MatrixXd chi0[3];
    chi0[0] = Eigen::MatrixXd(nel, 8);
    chi0[0].setZero();
    chi0[1] = Eigen::MatrixXd(nel, 8);
    chi0[1].setZero();
    chi0[2] = Eigen::MatrixXd(nel, 8);
    chi0[2].setZero();

    Eigen::MatrixXd chi0_e = Eigen::MatrixXd(8, 3);
    chi0_e.setZero();
    Eigen::MatrixXd ke = keMu + keLambda;
    Eigen::MatrixXd fe = feMu + feLambda;

    MatrixXd ke3_5end_3_5end(ke.rows() - 3, ke.cols() - 3);
    for (int i = 0; i < ke3_5end_3_5end.rows(); i++) {
        for (int j = 0; j < ke3_5end_3_5end.cols(); j++) {
            int indexi = i;
            int indexj = j;
            if (i == 0)
                indexi = 2;
            else
                indexi = i + 3;
            if (j == 0)
                indexj = 2;
            else
                indexj = j + 3;

            ke3_5end_3_5end(i, j) = ke(indexi, indexj);
        }
    }
    MatrixXd fe3_5end(fe.rows() - 3, fe.cols());
    for (int i = 0; i < fe3_5end.rows(); i++) {
        int indexi = i;
        if (i == 0)
            indexi = 2;
        else
            indexi = i + 3;

        fe3_5end.row(i) = fe.row(indexi);
    }

   /* std::cout << "ke3_5end_3_5end:\n" << ke3_5end_3_5end << "++++++++++++++++++++++\n";
    std::cout << "fe3_5end\n" << fe3_5end << "++++++++++++++++++++++\n";*/

    MatrixXd ke3_5cfe3_5 = ke3_5end_3_5end.inverse() * fe3_5end;

    for (int i = 0; i < ke3_5cfe3_5.rows(); i++) {
        int indexi = i;
        if (i == 0)
            indexi = 2;
        else
            indexi = i + 3;

        chi0_e.row(indexi) = ke3_5cfe3_5.row(i);
    }

    //std::cout << "chi0_e\n" << chi0_e << "++++++++++++++++++++++\n";

    VectorXd chi0_e_1 = chi0_e.col(0);
    VectorXd chi0_e_2 = chi0_e.col(1);
    VectorXd chi0_e_3 = chi0_e.col(2);

    VectorXd epsilonones = VectorXd::Ones(nel);

    chi0[0] = Eigen::kroneckerProduct(chi0_e_1.transpose(), epsilonones);
    chi0[1] = Eigen::kroneckerProduct(chi0_e_2.transpose(), epsilonones);
    chi0[2] = Eigen::kroneckerProduct(chi0_e_3.transpose(), epsilonones);

    //std::cout << "chi0[0]:\n" << chi0[0] << "++++++++++++++++++++++\n";
    //std::cout << "chi0[1]:\n" << chi0[1] << "++++++++++++++++++++++\n";
    //std::cout << "chi0[2]:\n" << chi0[2] << "++++++++++++++++++++++\n";

    //std::cout << "chi0_e1:\n" << chi0_e_1 << "++++++++++++++++++++++\n";
    //std::cout << "chi0_e2:\n" << chi0_e_2 << "++++++++++++++++++++++\n";
    //std::cout << "chi0_e3:\n" << chi0_e_3 << "++++++++++++++++++++++\n";

    Eigen::MatrixXd CH = Eigen::MatrixXd::Zero(3, 3);
    double cellVolume = lx * ly;
    Eigen::MatrixXd chiMatrix(chi);

    Eigen::MatrixXd chiMatrix2(chiMatrix.rows() + 2, chiMatrix.cols());
    chiMatrix2.row(0) = VectorXd::Zero(chiMatrix.cols());
    chiMatrix2.row(1) = VectorXd::Zero(chiMatrix.cols());
    chiMatrix2.block(2, 0, chiMatrix.rows(), chiMatrix.cols()) = chiMatrix;

    //std::cout << "chiMatrix:\n" << chiMatrix2 << "++++++++++++++++++++++\n";

    Map<VectorXd> chiVector(chiMatrix2.data(), chiMatrix2.size());

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Eigen::MatrixXd chiP1(edofMat.rows(), edofMat.cols());
            for (int indexx = 0; indexx < chiP1.rows(); indexx++) {
                for (int indexy = 0; indexy < chiP1.cols(); indexy++) {
                    chiP1(indexx, indexy) = chiVector(edofMat(indexx, indexy) + i * ndof - 1);
                }
            }
            Eigen::MatrixXd chiP2(edofMat.rows(), edofMat.cols());
            for (int indexx = 0; indexx < chiP1.rows(); indexx++) {
                for (int indexy = 0; indexy < chiP1.cols(); indexy++) {
                    chiP2(indexx, indexy) = chiVector(edofMat(indexx, indexy) + j * ndof - 1);
                }
            }

            //std::cout << "chiP1:\n" << chiP1 << "++++++++++++++++++++++\n";
            //std::cout << "Chi0:\n" << chi0[i] << "++++++++++++++++++++++\n";
            //std::cout << "keLambda:\n" << keLambda << "++++++++++++++++++++++\n";

            MatrixXd sumLambdaMult1 = (chi0[i] - chiP1) * keLambda;
            MatrixXd sumLambdaMult2 = chi0[j] - chiP2;
            //std::cout << "sumLambdaMult1:\n" << sumLambdaMult1 << "++++++++++++++++++++++\n";
            //std::cout << "sumLambdaMult2:\n" << sumLambdaMult2 << "++++++++++++++++++++++\n";

            MatrixXd sumLambdaT = ((sumLambdaMult1).cwiseProduct(sumLambdaMult2));
            MatrixXd sumMuT = (((chi0[i] - chiP1) * keMu).cwiseProduct(chi0[j] - chiP2));

            //std::cout << "sumLambdaT:\n" << sumLambdaT << "++++++++++++++++++++++\n";

            VectorXd sumLambdaRowsum = sumLambdaT.rowwise().sum();
            VectorXd sumMuRowsum = sumMuT.rowwise().sum();

            //std::cout << "sumLambdaRowWise:\n" << sumLambdaRowsum << "++++++++++++++++++++++\n";

            //VectorXd m;
            //m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
            //Matrix3d A = Map<Matrix3d>(m.data())

            MatrixXd sumLambda = Map<MatrixXd>(sumLambdaRowsum.data(), nely, nelx);
            MatrixXd sumMu = Map<MatrixXd>(sumMuRowsum.data(), nely, nelx);

            MatrixXd LambdaMatrix = Map<MatrixXd>(lambdaOfx.data(), nely, nelx);
            MatrixXd MuMatrix = Map<MatrixXd>(muOfx.data(), nely, nelx);

            /*std::cout << "iter  "<<i*3+j<<endl<<"sumLambda:\n" << sumLambda << "++++++++++++++++++++++\n";
            std::cout << "sumMu:\n" << sumMu << "++++++++++++++++++++++\n";
            std::cout << "LambdaMatrix:\n" << LambdaMatrix << "++++++++++++++++++++++\n";
            std::cout << LambdaMatrix.cwiseProduct(sumLambda) << "++++++++++++++++++++++\n";
            std::cout << cellVolume<<endl<<(LambdaMatrix.cwiseProduct(sumLambda) + MuMatrix.cwiseProduct(sumMu)).sum() << "++++++++++++++++++++++\n";*/

            //std::cout << "LambdaVec.*sumLambda:\n" << LambdaMatrix.cwiseProduct(sumLambda) << "+++++++++++++++++++\n";
            //std::cout << "MuVec.*sumMu:\n" << MuMatrix.cwiseProduct(sumMu) << "+++++++++++++++++++\n";

            CH(i, j) = 1.0 / cellVolume
                * (LambdaMatrix.cwiseProduct(sumLambda) + MuMatrix.cwiseProduct(sumMu)).sum();
            //cout << "CH:"<<endl<<CH << endl;
        }
    }

    //std::cout << "CH:\n" << CH << "++++++++++++++++++++++\n";

    return CH;
}

void elementMatVec(double a,
    double b,
    double phi,
    MatrixXd& keLambda,
    MatrixXd& keMu,
    MatrixXd& feLambda,
    MatrixXd& feMu)
{
    // Constitutive matrix contributions
    Matrix3d CMu = Matrix3d::Zero();
    CMu.diagonal() << 2, 2, 1;
    Matrix3d CLambda = Matrix3d::Zero();
    for (int CLindexx = 0; CLindexx < 2; CLindexx++) {
        for (int CLindexy = 0; CLindexy < 2; CLindexy++) {
            CLambda(CLindexx, CLindexy) = 1;
        }
    }

    //std::cout << "CLambda:::\n" << CLambda << '\n';

    // Two Gauss points in both directions
    double xx[2] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    double yy[2] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    double ww[2] = { 1.0, 1.0 };

    // Initialize
    keLambda = MatrixXd::Zero(8, 8);
    keMu = MatrixXd::Zero(8, 8);
    feLambda = MatrixXd::Zero(8, 3);
    feMu = MatrixXd::Zero(8, 3);

    Matrix<double, 3, 4> L;
    L << 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0;

    for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 2; jj++) {
            // Integration point
            double x = xx[ii];
            double y = yy[jj];

            // Differentiated shape functions
            VectorXd dNx(4), dNy(4);
            dNx << -(1 - y), 1 - y, 1 + y, -(1 + y);
            dNy << -(1 - x), -(1 + x), 1 + x, 1 - x;
            dNx /= 4.0;
            dNy /= 4.0;

            //std::cout << dNx << '\n';
            //std::cout << dNy << '\n';

            // Jacobian
            MatrixXd dN_t(2, 4);
            dN_t.row(0) = dNx;
            dN_t.row(1) = dNy;

            MatrixXd dN_tt(4, 2);
            VectorXd vec1(4), vec2(4);
            vec1 << -a, a, a + 2 * b / tan(phi * M_PI / 180), 2 * b / tan(phi * M_PI / 180) - a;
            vec2 << -b, -b, b, b;
            dN_tt.col(0) = vec1;
            dN_tt.col(1) = vec2;

            //std::cout << "dN_t:::\n" << dN_t << '\n';
            //std::cout << "dN_tt:::\n" << dN_tt << '\n';

            Matrix2d J = dN_t * dN_tt;
            //J.row(0) << dN_t.row(0) * Matrix<double, 4, 1>(-a, a, a + 2 * b / tan(phi * M_PI / 180), 2 * b / tan(phi * M_PI / 180) - a);
            //J.row(1) << dN_t.row(1) * Matrix<double, 4, 1>(-b, -b, b, b);

            //std::cout << "J:::\n" << J << '\n';

            double detJ = J.determinant();
            Matrix2d invJ;
            invJ << J(1, 1), -J(0, 1), -J(1, 0), J(0, 0);
            invJ *= 1.0 / detJ;

            // Weight factor at this point
            double weight = ww[ii] * ww[jj] * detJ;

            // Strain-displacement matrix
            MatrixXd G(4, 4);
            G.setZero();
            G.block<2, 2>(0, 0) = invJ;
            G.block<2, 2>(2, 2) = invJ;
            MatrixXd dN = Matrix<double, 4, 8>::Zero();
            for (int indexOfdN = 0; indexOfdN < 4; indexOfdN++) {
                dN(0, indexOfdN * 2) = dNx(indexOfdN);
                dN(1, indexOfdN * 2) = dNy(indexOfdN);
                dN(2, indexOfdN * 2 + 1) = dNx(indexOfdN);
                dN(3, indexOfdN * 2 + 1) = dNy(indexOfdN);
            }

            MatrixXd B = L * G * dN;
            //std::cout << "G:::\n" << G << '\n';
            //std::cout << "dN:::\n" << dN << '\n';
            //std::cout << "B:::\n" << B << '\n';

            // Element matrices
            keLambda += weight * (B.transpose() * CLambda * B);
            keMu += weight * (B.transpose() * CMu * B);
            // Element loads
            feLambda += weight * (B.transpose() * CLambda * MatrixXd::Identity(3, 3));
            feMu += weight * (B.transpose() * CMu * MatrixXd::Identity(3, 3));

            //std::cout << keLambda << '\n';
        }
    }
}


MatrixXd homogenize_therm(
    int lx, int ly, std::vector<double> lambda, std::vector<double> mu, double phi, MatrixXi x)
{
    //INITIALIZE
    // Deduce discretization
    int nely = x.rows(), nelx = x.cols();
    // Stiffness matrix consists of two parts, one belonging to lambda and
    // one belonging to mu. Same goes for load vector
    double dx = (double)lx / nelx, dy = (double)ly / nely;
    int nel = nelx * nely;
    MatrixXd keLambda, keMu, feLambda, feMu;
    elementMatVec_therm(dx / 2, dy / 2, phi, keLambda, keMu, feLambda, feMu);
    for (int i = 0; i < keMu.rows(); i += 2) {      // 遍历奇数行（MATLAB索引1:2:end对应C++索引0,2,4...）
        for (int j = 0; j < keMu.cols(); j += 2) {  // 遍历奇数列
            if (i + 1 < keMu.rows() && j + 1 < keMu.cols()) { // 防止越界访问
                keMu(i, j) += keMu(i + 1, j + 1);       // 偶位置元素累加到奇位置
            }
        }
    }
    //cout << keMu << endl<< "kela:"<<keLambda<<"fela:"<<feLambda<<"fem:"<<feMu<<endl;
    // Node numbers and element degrees of freedom for full (not periodic) mesh
    //nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
    MatrixXi nodenrs(nely + 1, nelx + 1);
    for (int i = 0, c = 1; i <= nelx; ++i) {
        for (int j = 0; j <= nely; ++j) {
            nodenrs(j, i) = c++;
        }
    }
    //cout << nodenrs.rows() << "   " << nodenrs.cols() << endl << nodenrs << endl;
    //edofVec = reshape(2 * nodenrs(1:end - 1, 1 : end - 1) + 1, nel, 1);
    VectorXi edofVec(nel); 
     
    for (int i = 0, c = 0; i < nelx; ++i) {
        for (int j = 0; j < nely; ++j) {
            edofVec(c) = 2 * nodenrs(j, i) + 1;
            c++;
        }
    }
    //cout << edofVec << endl;
    //edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2 * nely + [2 3 0 1] - 2 - 1], nel, 1);
    MatrixXi edofMat(nel, 8);
    edofMat.col(0) = edofVec;
    edofMat.col(1) = edofVec.array() + 1;
    edofMat.col(2) = edofVec.array() + 2 * nely + 2;
    edofMat.col(3) = edofVec.array() + 2 * nely + 3;
    edofMat.col(4) = edofVec.array() + 2 * nely;
    edofMat.col(5) = edofVec.array() + 2 * nely + 1;
    edofMat.col(6) = edofVec.array() - 2;
    edofMat.col(7) = edofVec.array() - 1;
  
    //--------------------------Test----------------------------------------
    //std::cout << "keLambda:\n" << keLambda << "++++++++++++++++++++++\n";
    //std::cout << "keMu:\n" << keMu << "++++++++++++++++++++++\n";
    //std::cout << "feLambda:\n" << feLambda << "++++++++++++++++++++++\n";
    //std::cout << "feMu:\n" << feMu << "++++++++++++++++++++++\n";

    //std::cout << "nodenrs:\n" << nodenrs << "++++++++++++++++++++++\n";
    //std::cout << "edofVec:\n" << edofVec << "++++++++++++++++++++++\n";
    //std::cout << "edofMat:\n" << edofMat << "++++++++++++++++++++++\n";
    //----------------------------------------------------------------------
    //IMPOSE PERIODIC BOUNDARY CONDITIONS
    //Use original edofMat to index into list with the periodic dofs
    int nn = (nelx + 1) * (nely + 1); // Total number of nodes
    int nnP = nelx * nely;            // Total number of unique nodes
    Eigen::MatrixXi nnPArray(nely + 1, nelx + 1);
    nnPArray.setZero();
    for (int i = 0; i < nely; i++)
        for (int j = 0; j < nelx; j++)
            nnPArray(j, i) = i * nelx + j + 1;
    // Extend with a mirror of the top border
    for (int j = 0; j < nelx + 1; j++) {
        nnPArray(nely, j) = nnPArray(0, j);
    }
    // Extend with a mirror of the left border
    for (int i = 0; i < nely + 1; i++) {
        nnPArray(i, nelx) = nnPArray(i, 0);
    }
    Eigen::VectorXd nnParrayVector(nn);
    int indexOfnnParrayVector = 0;
    for (int i = 0; i <= nelx; i++) {
        for (int j = 0; j <= nely; j++) {
            nnParrayVector(indexOfnnParrayVector++) = nnPArray(j, i);
        }
    }
    //cout <<"nnParrayVector:"<< nnParrayVector << endl;
    Eigen::VectorXd dofVector(2 * nn);
    dofVector.setZero();

    for (int i = 0; i < nn; i++) {
        dofVector(i * 2) = 2 * nnParrayVector(i) - 1;
        dofVector(i * 2 + 1) = 2 * nnParrayVector(i);
    }

    //for (int i = 0; i <= nely; i++)
    //{
    //	for (int j = 0; j <= nelx; j++)
    //	{
    //		int p = i * nelx + nelx - j;
    //		dofVector(2 * p) = 2 * (nnPArray(j, i)) - 1;
    //		dofVector(2 * p + 1) = 2 * (nnPArray(j, i));
    //	}
    //}

    Eigen::MatrixXi edofMatNew(nel, 8);
    for (int i = 0; i < nel; i++)
        for (int j = 0; j < 8; j++)
        {
            edofMatNew(i, j) = dofVector(edofMat(i, j) - 1);
            //cout<<i<<"  "<<j<<"   edofMat(i, j) - 1:"<< edofMat(i, j) - 1<<"     val:"<< dofVector(edofMat(i, j) - 1) << endl;
        }
    edofMat = edofMatNew;
    int ndof = 2 * nnP; // Number of dofs
    //std::cout << "edofMat:\n" << ndof<<endl<<edofMat.cols() << endl << edofMat.row(0) << "++++++++++++++++++++++\n";
    //--------------------------Test----------------------------------------
    //std::cout << "nnPArray:\n" << nnPArray << "++++++++++++++++++++++\n";
    //std::cout << "nnParrayVector:\n" << nnParrayVector << "++++++++++++++++++++++\n";
    //std::cout << "dofVector:\n" << dofVector << "++++++++++++++++++++++\n";
    //std::cout << "edofMat:\n" << edofMat << "++++++++++++++++++++++\n";
    //std::cout << "feMu:\n" << feMu << "++++++++++++++++++++++\n";
    //std::cout << "nodenrs:\n" << nodenrs << "++++++++++++++++++++++\n";
    //std::cout << "edofVec:\n" << edofVec << "++++++++++++++++++++++\n";
    //std::cout << "edofMat:\n" << edofMat << "++++++++++++++++++++++\n";
    //----------------------------------------------------------------------

    // ASSEMBLE STIFFNESS MATRIX
    // Indexing vectors
    //Eigen::MatrixXi iK = Eigen::Map<Eigen::MatrixXi>(edofMat.data(), edofMat.size(), 1).replicate(1, 8).transpose();
    //Eigen::MatrixXi jK = Eigen::Map<Eigen::MatrixXi>(edofMat.data(), edofMat.size(), 1).replicate(8, 1);

    Eigen::MatrixXi iKones = Eigen::MatrixXi::Ones(8, 1);
    Eigen::MatrixXi iK = Eigen::kroneckerProduct(edofMat, iKones).transpose();

    Eigen::MatrixXi jKones = Eigen::MatrixXi::Ones(1, 8);
    Eigen::MatrixXi jK = Eigen::kroneckerProduct(edofMat, jKones).transpose();

    //std::cout << "iK:\n" << iK.size()<<endl<<iK.col(0)<< "++++++++++++++++++++++\n";
    //std::cout << "jK:\n" << jK.col(0) << "++++++++++++++++++++++\n";

    // Material properties in the different elements
    Eigen::VectorXd lambdaOfx(nel);
    Eigen::VectorXd muOfx(nel);
    for (int i = 0; i < nel; i++) {
        if (x(i) == 1) {
            lambdaOfx(i) = lambda[0];
            muOfx(i) = mu[0];
        }
        else if (x(i) == 2) {
            lambdaOfx(i) = lambda[1];
            muOfx(i) = mu[1];
        }
    }
    // The corresponding stiffness matrix entries
    Eigen::VectorXd keLambdaVector(keLambda.rows() * keLambda.cols());
    Eigen::VectorXd keMuVector(keMu.rows() * keMu.cols());
    int indexOfkeLambdaVector = 0;
    for (int i = 0; i < keLambda.cols(); i++) {
        for (int j = 0; j < keLambda.rows(); j++) {
            keLambdaVector(indexOfkeLambdaVector++) = keLambda(j, i);
        }
    }
    int indexOfkeMuVector = 0;
    for (int i = 0; i < keMu.cols(); i++) {
        for (int j = 0; j < keMu.rows(); j++) {
            keMuVector(indexOfkeMuVector++) = keMu(j, i);
        }
    }
    Eigen::MatrixXd sK = (keLambdaVector * lambdaOfx.transpose() + keMuVector * muOfx.transpose());

    //std::cout << "sK:\n" << sK << "++++++++++++++++++++++\n";

    //std::cout << iK.rows() * iK.cols() << "   ===    " << jK.rows() * jK.cols() << "  ===  " << sK.rows() * sK.cols() << '\n';
    //Map<RowVectorXf> v1(M1.data(), M1.size());
    //cout << "v1:" << endl << v1 << endl;

    //Matrix<float, Dynamic, Dynamic, RowMajor> M2(M1);
    //Map<RowVectorXf> v2(M2.data(), M2.size());
    //cout << "v2:" << endl << v2 << endl;

    Map<VectorXi> iKVector(iK.data(), iK.size());
    Map<VectorXi> jKVector(jK.data(), jK.size());
    Map<VectorXd> sKVector(sK.data(), sK.size());

    Eigen::SparseMatrix<double> K = Eigen::SparseMatrix<double>(ndof, ndof);
    //// Fill K with triplets
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < iK.size(); ++i) {
        triplets.emplace_back(iKVector(i) - 1, jKVector(i) - 1, sKVector(i));
    }
    K.setFromTriplets(triplets.begin(), triplets.end());

    //--------------------------Test----------------------------------------
    //----------------------------------------------------------------------

    //std::cout << "K:\n" << MatrixXd(K) << "++++++++++++++++++++++\n";

    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //LOAD VECTORS AND SOLUTION
    Eigen::VectorXd feLambdaVector(feLambda.rows() * feLambda.cols());
    Eigen::VectorXd feMuVector(feMu.rows() * feMu.cols());
    int indexOffeLambdaVector = 0;
    for (int i = 0; i < feLambda.cols(); i++) {
        for (int j = 0; j < feLambda.rows(); j++) {
            feLambdaVector(indexOffeLambdaVector++) = feLambda(j, i);
        }
    }
    int indexOffeMuVector = 0;
    for (int i = 0; i < feMu.cols(); i++) {
        for (int j = 0; j < feMu.rows(); j++) {
            feMuVector(indexOffeMuVector++) = feMu(j, i);
        }
    }

    Eigen::MatrixXd sF = (feLambdaVector * lambdaOfx.transpose() + feMuVector * muOfx.transpose());
    //std::cout << "feLambdaVector:\n" << feLambdaVector << "++++++++++++++++++++++\n";
    //std::cout << "feMuVector:\n" << feMuVector << "++++++++++++++++++++++\n";
    //std::cout << "sF:\n" << sF << "++++++++++++++++++++++\n";

    //			//nel * 8
    Eigen::MatrixXi iF(24, nel);
    Eigen::MatrixXi edofMatT = edofMat.transpose();
    iF.block(0, 0, edofMatT.rows(), edofMatT.cols()) = edofMatT;
    iF.block(edofMatT.rows(), 0, edofMatT.rows(), edofMatT.cols()) = edofMatT;
    iF.block(edofMatT.rows() * 2, 0, edofMatT.rows(), edofMatT.cols()) = edofMatT;

    Eigen::MatrixXi jF(24, nel);
    Eigen::MatrixXi jFzeros = Eigen::MatrixXi::Ones(8, nel);
    jF.block(0, 0, jFzeros.rows(), jFzeros.cols()) = jFzeros;
    jF.block(jFzeros.rows(), 0, jFzeros.rows(), jFzeros.cols()) = jFzeros * 2;
    jF.block(jFzeros.rows() * 2, 0, jFzeros.rows(), jFzeros.cols()) = jFzeros * 3;

    Map<VectorXi> iFVector(iF.data(), iF.size());
    Map<VectorXi> jFVector(jF.data(), jF.size());
    Map<VectorXd> sFVector(sF.data(), sF.size());

    Eigen::SparseMatrix<double> F(ndof, 3);

    std::vector<Eigen::Triplet<double>> triplets2;
    for (int i = 0; i < iF.size(); ++i) {
        triplets2.emplace_back(iFVector(i) - 1, jFVector(i) - 1, sFVector(i));
    }
    F.setFromTriplets(triplets2.begin(), triplets2.end());
    //showSparseMatrix(F.block(0,0,100,3));
    //chi(3:2:ndof,:) = K(3:2:end,3:2:end)\[F(3:2:end,1) F(4:2:end,2)];
    // 生成子矩阵的行索引 (MATLAB: 3:2:end → C++索引 2,4,6...)
    std::vector<int> sub_rows;
    for (int i = 2; i < ndof; i += 2) {
        sub_rows.push_back(i);
    }
    // 提取 K 的数据
    int sub_size = sub_rows.size();
    Eigen::SparseMatrix<double> K_3_2_end(sub_size, sub_size); // 子矩阵大小
    std::vector<Eigen::Triplet<double>> K_sub_triplets;

    // 遍历 K 中所有非零元素，筛选满足行列都在 sub_indices 中的元素
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            int row = it.row();
            int col = it.col();

            // 判断 row 和 col 是否都在 sub_indices 中
            auto it_row = std::find(sub_rows.begin(), sub_rows.end(), row);
            auto it_col = std::find(sub_rows.begin(), sub_rows.end(), col);

            if (it_row != sub_rows.end() && it_col != sub_rows.end()) {
                int new_row = std::distance(sub_rows.begin(), it_row);
                int new_col = std::distance(sub_rows.begin(), it_col);
                K_sub_triplets.emplace_back(new_row, new_col, it.value());
            }
        }
    }
    K_3_2_end.setFromTriplets(K_sub_triplets.begin(), K_sub_triplets.end());
    
    // 提取 F(3:2:end, 1) => F(2,0), F(4,0), ...  提取 F(4:2:end, 2) => F(3,1), F(5,1), ...
    std::vector<Eigen::Triplet<double>> F_sub_triplets;
    int row_sub = 0;

    // 第1列：F(3:2:end,1) ⇒ 索引 2,4,6,...
    for (int i = 2; i < ndof; i += 2) {
        double val = F.coeff(i, 0);
        if (val != 0.0) {
            F_sub_triplets.emplace_back(row_sub, 0, val);
        }
        ++row_sub;
    }
    int len1 = row_sub;
    row_sub = 0;
    // 第2列：F(4:2:end,2) ⇒ 索引 3,5,7,...
    for (int i = 3; i < ndof; i += 2) {
        double val = F.coeff(i, 1);
        if (val != 0.0) {
            F_sub_triplets.emplace_back(row_sub, 1, val);
        }
        ++row_sub;
    }
    assert(len1 == row_sub);

    // 构建最终稀疏矩阵 F_sub（max_row x 2）
    Eigen::SparseMatrix<double> F_sub(len1, 2);
    F_sub.setFromTriplets(F_sub_triplets.begin(), F_sub_triplets.end());
    
    SparseLU<SparseMatrix<double>> solver;
    //solver.compute(K.transpose().block(2, 2,  ndof - 2, ndof - 2)); // 求解K中3到ndof的部分，注意要去除边界上的自由度
    solver.compute(K_3_2_end); // 求解K_3_2_end
    
    //chi = solver.solve(F_sub); // 求解线性方程组
    //cout << "chi computed!" << endl;
    Eigen::SparseMatrix<double>  X1 = solver.solve(F_sub); // 求解线性方程组
    //showSparseMatrix(X1);

    // 将解结果填充到 chi 的稀疏矩阵中（仅非零值）
    std::vector<int> target_rows;
    for (int i = 2; i < ndof; i += 2) { // MATLAB 3:2:ndof == C++ 2,4,6,...
        target_rows.push_back(i);
    }

    Eigen::SparseMatrix<double> chi(ndof, 2);
    int num_target_rows = target_rows.size();
    // 将 X1 的每个非零元素赋值到 chi 的对应奇数行（从第3行起）
    for (int k = 0; k < X1.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(X1, k); it; ++it) {
            int row_in_X1 = it.row(); // 行号 in X1
            int col = it.col();
            int row_in_chi = target_rows[row_in_X1];
            chi.coeffRef(row_in_chi, col) = it.value();
        }
    }
    cout << "chi computed!" <<endl;
    //showSparseMatrix(chi);

//HOMOGENIZATION
//注意问题  matlab的下标可能与Eigen不同

    Eigen::MatrixXd chi0[3];
    chi0[0] = Eigen::MatrixXd(nel, 8);
    chi0[0].setZero();
    chi0[1] = Eigen::MatrixXd(nel, 8);
    chi0[1].setZero();
    chi0[2] = Eigen::MatrixXd(nel, 8);
    chi0[2].setZero();

    Eigen::MatrixXd chi0_e = Eigen::MatrixXd(8, 3);
    chi0_e.setZero();
   /* Eigen::MatrixXd ke = keMu + keLambda;
    Eigen::MatrixXd fe = feMu + feLambda;*/

    //chi0_e(3:2:end,1:2) = keMu(3:2:end,3:2:end)\[feMu(3:2:end,1) feMu(4:2:end,2)];
    std::vector<int> kemu_rows_cols; // kemu目标行列索引 
    std::vector<int> kemu_rows_cols2; // kemu目标行列索引 
    for (int i = 2; i < keMu.rows(); i += 2) {
        kemu_rows_cols.push_back(i);
    }

    for (int i = 2; i < keMu.cols(); i += 2) {
        kemu_rows_cols2.push_back(i);
    }
    //keMu(3:2:end,3:2:end)
    MatrixXd keMu3_2_end_3_2_end(kemu_rows_cols.size(), kemu_rows_cols.size());
    for (int i = 0; i < keMu3_2_end_3_2_end.rows(); i++) {
        for (int j = 0; j < keMu3_2_end_3_2_end.cols(); j++) {
            keMu3_2_end_3_2_end(i, j) = keMu(kemu_rows_cols[i], kemu_rows_cols[j]);
        }
    }
   // cout << "keMu3_2_end_3_2_end: " << keMu3_2_end_3_2_end << endl;
    //[F(3:2:end,1) F(4:2:end,2)]
    std::vector<int> femu_rows1; // feMu奇数行第一列行列索引 
    std::vector<int> femu_rows2; // feMu偶数行第二列行列索引 
    for (int i = 2; i < feMu.rows(); i += 2) {
        femu_rows1.push_back(i);
    }
    for (int i = 3; i < feMu.rows(); i += 2) {
        femu_rows2.push_back(i);
    }
    assert(femu_rows1.size() == femu_rows2.size());

    MatrixXd Femu_sub(femu_rows1.size(), 2);
    for (int i = 0; i < femu_rows1.size(); i++) {
        Femu_sub(i, 0) = feMu(femu_rows1[i], 0); // 第一列
        Femu_sub(i, 1) = feMu(femu_rows2[i], 1); // 第二列
    }
    // chi0_e(3:2:end,1:2)
    Eigen::FullPivLU<Eigen::MatrixXd> lu_solver;
    lu_solver.compute(keMu3_2_end_3_2_end);
    Eigen::MatrixXd X2 = lu_solver.solve(Femu_sub);
    cout << X2 << endl;
    // ===================================================================
    for (int i = 0; i < kemu_rows_cols.size(); ++i) {
        const int row = kemu_rows_cols[i];
        chi0_e(row, 0) = X2(i, 0);  // 第一列
        chi0_e(row, 1) = X2(i, 1);  // 第二列
    }

    //std::cout << "chi0_e\n" << chi0_e << "++++++++++++++++++++++\n";

    VectorXd chi0_e_1 = chi0_e.col(0);
    VectorXd chi0_e_2 = chi0_e.col(1);
    VectorXd chi0_e_3 = chi0_e.col(2);

    VectorXd epsilonones = VectorXd::Ones(nel);

    chi0[0] = Eigen::kroneckerProduct(chi0_e_1.transpose(), epsilonones);
    chi0[1] = Eigen::kroneckerProduct(chi0_e_2.transpose(), epsilonones);
    chi0[2] = Eigen::kroneckerProduct(chi0_e_3.transpose(), epsilonones);

    //std::cout << "chi0[0]:\n" << chi0[0] << "++++++++++++++++++++++\n";
    //std::cout << "chi0[1]:\n" << chi0[1] << "++++++++++++++++++++++\n";
    //std::cout << "chi0[2]:\n" << chi0[2] << "++++++++++++++++++++++\n";

    //std::cout << "chi0_e1:\n" << chi0_e_1 << "++++++++++++++++++++++\n";
    //std::cout << "chi0_e2:\n" << chi0_e_2 << "++++++++++++++++++++++\n";
    //std::cout << "chi0_e3:\n" << chi0_e_3 << "++++++++++++++++++++++\n";

    Eigen::MatrixXd CH = Eigen::MatrixXd::Zero(3, 3);
    double cellVolume = lx * ly;
    Eigen::MatrixXd chiMatrix(chi);
    Map<VectorXd> chiVector(chiMatrix.data(), chiMatrix.size());

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            Eigen::MatrixXd chiP1(edofMat.rows(), edofMat.cols());
            Eigen::MatrixXd chiP2(edofMat.rows(), edofMat.cols());
            for (int indexx = 0; indexx < chiP1.rows(); indexx++) {
                for (int indexy = 0; indexy < chiP1.cols(); indexy++) {
                    chiP1(indexx, indexy) = chiVector(edofMat(indexx, indexy) + i * ndof - 1);
                    chiP2(indexx, indexy) = chiVector(edofMat(indexx, indexy) + j * ndof - 1);
                }
            }
            
           /* for (int indexx = 0; indexx < chiP1.rows(); indexx++) {
                for (int indexy = 0; indexy < chiP1.cols(); indexy++) {
                    chiP2(indexx, indexy) = chiVector(edofMat(indexx, indexy) + j * ndof - 1);
                }
            }*/

            //std::cout << "chiP1:\n" << chiP1 << "++++++++++++++++++++++\n";
            //std::cout << "Chi0:\n" << chi0[i] << "++++++++++++++++++++++\n";
            //std::cout << "keLambda:\n" << keLambda << "++++++++++++++++++++++\n";
            
            MatrixXd sumLambdaMult1 = (chi0[i] - chiP1) * keLambda;
            MatrixXd sumLambdaMult2 = chi0[j] - chiP2;
            //std::cout << "sumLambdaMult1:\n" << sumLambdaMult1 << "++++++++++++++++++++++\n";
            //std::cout << "sumLambdaMult2:\n" << sumLambdaMult2 << "++++++++++++++++++++++\n";

            MatrixXd sumLambdaT = ((sumLambdaMult1).cwiseProduct(sumLambdaMult2));
            MatrixXd sumMuT = (((chi0[i] - chiP1) * keMu).cwiseProduct(chi0[j] - chiP2));
            //std::cout << "sumLambdaT:\n" << sumLambdaT << "++++++++++++++++++++++\n";

            VectorXd sumLambdaRowsum = sumLambdaT.rowwise().sum();
            VectorXd sumMuRowsum = sumMuT.rowwise().sum();


            MatrixXd sumLambda = Map<MatrixXd>(sumLambdaRowsum.data(), nely, nelx);
            MatrixXd sumMu = Map<MatrixXd>(sumMuRowsum.data(), nely, nelx);

            MatrixXd LambdaMatrix = Map<MatrixXd>(lambdaOfx.data(), nely, nelx);
            MatrixXd MuMatrix = Map<MatrixXd>(muOfx.data(), nely, nelx);

            //std::cout << "sumLambda:\n" << sumLambda << "++++++++++++++++++++++\n";
            //std::cout << "sumMu:\n" << sumMu << "++++++++++++++++++++++\n";
            /*std::cout << "LambdaMatrix:\n" << LambdaMatrix << "++++++++++++++++++++++\n";
            std::cout << "MuMatrix:\n" << MuMatrix << "++++++++++++++++++++++\n";*/

            //std::cout << "LambdaVec.*sumLambda:\n" << LambdaMatrix.cwiseProduct(sumLambda) << "+++++++++++++++++++\n";
            //std::cout << "MuVec.*sumMu:\n" << MuMatrix.cwiseProduct(sumMu) << "+++++++++++++++++++\n";
            CH(i, j) = 1.0 / cellVolume
                * (LambdaMatrix.cwiseProduct(sumLambda) + MuMatrix.cwiseProduct(sumMu)).sum();
        }
    }
    
    std::cout << "CH:\n" << CH <<endl<< "++++++++++++++++++++++\n";
    Eigen::MatrixXd Kappa_H = CH.block(0, 0, 2, 2);
    return Kappa_H;
}

void elementMatVec_therm(double a,
    double b,
    double phi,
    MatrixXd& keLambda,
    MatrixXd& keMu,
    MatrixXd& feLambda,
    MatrixXd& feMu)
{
    // Constitutive matrix contributions
    Matrix3d CMu = Matrix3d::Zero();
    CMu.diagonal() << 1, 1, 0;
    Matrix3d CLambda = Matrix3d::Zero();


    // Two Gauss points in both directions
    double xx[2] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    double yy[2] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    double ww[2] = { 1.0, 1.0 };

    // Initialize
    keLambda = MatrixXd::Zero(8, 8);
    keMu = MatrixXd::Zero(8, 8);
    feLambda = MatrixXd::Zero(8, 3);
    feMu = MatrixXd::Zero(8, 3);

    Matrix<double, 3, 4> L;
    L << 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0;

    for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 2; jj++) {
            // Integration point
            double x = xx[ii];
            double y = yy[jj];

            // Differentiated shape functions
            VectorXd dNx(4), dNy(4);
            dNx << -(1 - y), 1 - y, 1 + y, -(1 + y);
            dNy << -(1 - x), -(1 + x), 1 + x, 1 - x;
            dNx /= 4.0;
            dNy /= 4.0;

            //std::cout << dNx << '\n';
            //std::cout << dNy << '\n';

            // Jacobian
            MatrixXd dN_t(2, 4);
            dN_t.row(0) = dNx;
            dN_t.row(1) = dNy;

            MatrixXd dN_tt(4, 2);
            VectorXd vec1(4), vec2(4);
            vec1 << -a, a, a + 2 * b / tan(phi * M_PI / 180), 2 * b / tan(phi * M_PI / 180) - a;
            vec2 << -b, -b, b, b;
            dN_tt.col(0) = vec1;
            dN_tt.col(1) = vec2;

            //std::cout << "dN_t:::\n" << dN_t << '\n';
            //std::cout << "dN_tt:::\n" << dN_tt << '\n';

            Matrix2d J = dN_t * dN_tt;
            //J.row(0) << dN_t.row(0) * Matrix<double, 4, 1>(-a, a, a + 2 * b / tan(phi * M_PI / 180), 2 * b / tan(phi * M_PI / 180) - a);
            //J.row(1) << dN_t.row(1) * Matrix<double, 4, 1>(-b, -b, b, b);

            //std::cout << "J:::\n" << J << '\n';

            double detJ = J.determinant();
            Matrix2d invJ;
            invJ << J(1, 1), -J(0, 1), -J(1, 0), J(0, 0);
            invJ *= 1.0 / detJ;

            // Weight factor at this point
            double weight = ww[ii] * ww[jj] * detJ;

            // Strain-displacement matrix
            MatrixXd G(4, 4);
            G.setZero();
            G.block<2, 2>(0, 0) = invJ;
            G.block<2, 2>(2, 2) = invJ;
            MatrixXd dN = Matrix<double, 4, 8>::Zero();
            for (int indexOfdN = 0; indexOfdN < 4; indexOfdN++) {
                dN(0, indexOfdN * 2) = dNx(indexOfdN);
                dN(1, indexOfdN * 2) = dNy(indexOfdN);
                dN(2, indexOfdN * 2 + 1) = dNx(indexOfdN);
                dN(3, indexOfdN * 2 + 1) = dNy(indexOfdN);
            }

            MatrixXd B = L * G * dN;
            //std::cout << "G:::\n" << G << '\n';
            //std::cout << "dN:::\n" << dN << '\n';
            //std::cout << "B:::\n" << B << '\n';

            // Element matrices
            keLambda += weight * (B.transpose() * CLambda * B);
            keMu += weight * (B.transpose() * CMu * B);
            // Element loads
            feLambda += weight * (B.transpose() * CLambda * MatrixXd::Identity(3, 3));
            feMu += weight * (B.transpose() * CMu * MatrixXd::Identity(3, 3));

            //std::cout << keLambda << '\n';
        }
    }
}


void showSparseMatrix(SparseMatrix<double> X)
{
    cout << X.rows() << "x" << X.cols() << endl;
    for (int k = 0; k < X.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(X, k); it; ++it) {
            std::cout << "(" << it.row() + 1 << ", " << it.col() + 1 << ") "
                << it.value() << std::endl;
        }
    }
}


//math tools
double BSampleFunction::F03(double t)
{
    return 1.0 / 6 * (-t * t * t + 3 * t * t - 3 * t + 1);
}
double BSampleFunction::F13(double t)
{
    return 1.0 / 6 * (3 * t * t * t - 6 * t * t + 4);
}
double BSampleFunction::F23(double t)
{
    return 1.0 / 6 * (-3 * t * t * t + 3 * t * t + 3 * t + 1);
}
double BSampleFunction::F33(double t)
{
    return 1.0 / 6 * t * t * t;
}

std::vector<Eigen::Vector2d> BSampleFunction::ThreeOrderBSplineInterpolatePt(std::vector<Eigen::Vector2d>& pt, int InsertNum)
{
    if (pt.size() == 0 || InsertNum <= 0)
        return std::vector<Eigen::Vector2d>();
    int Num = pt.size();
    int InsertNumSum = 0;
    for (int i = 0; i < Num - 1; i++)
        InsertNumSum += InsertNum;

    std::vector<Eigen::Vector2d> temp(Num + 2);
    for (int i = 0; i < Num; i++)
        temp[i + 1] = pt[i];

    temp[0] = Eigen::Vector2d(2 * temp[1].x() - temp[2].x(), 2 * temp[1].y() - temp[2].y());
    temp[Num + 1] = Eigen::Vector2d(2 * temp[Num].x() - temp[Num - 1].x(),
        2 * temp[Num].y() - temp[Num - 1].y());


    Eigen::Vector2d NodePt1, NodePt2, NodePt3, NodePt4;
    double t;
    std::vector<Eigen::Vector2d> newpt(Num + InsertNumSum);

    int totalnum = 0;
    for (int i = 0; i < Num - 1; i++)
    {
        NodePt1 = temp[i];
        NodePt2 = temp[i + 1];
        NodePt3 = temp[i + 2];
        NodePt4 = temp[i + 3];
        double dt = 1.0 / (InsertNum + 1);

        for (int j = 0; j < InsertNum + 1; j++) {
            t = dt * j;
            newpt[totalnum] = F03(t) * NodePt1 + F13(t) * NodePt2 + F23(t) * NodePt3
                + F33(t) * NodePt4;
            totalnum++;
        }

        if (i == Num - 2) {
            t = 1;
            newpt[totalnum] = F03(t) * NodePt1 + F13(t) * NodePt2 + F23(t) * NodePt3
                + F33(t) * NodePt4;
            totalnum++;
        }
    }

    return newpt;

    Num = Num + InsertNumSum;
}

HermitFunction::HermitFunction(const Eigen::MatrixXd& control_points, const Eigen::MatrixXd& tangents)
    : m_P(control_points)
    , m_T(tangents)
{
    assert(control_points.cols() == tangents.cols() && control_points.rows() == tangents.rows());
}

HermitFunction::HermitFunction(Eigen::Vector2d b, Eigen::Vector2d p0, Eigen::Vector2d e, Eigen::Vector2d p1)
{
    Eigen::MatrixXd P(2, 2);
    P << b.x(), e.x(), b.y(), e.y();
    m_P = P;

    Eigen::MatrixXd T(2, 2);
    T << p0.x(), p1.x(), p0.y(), p1.y();
    m_T = T;
}

HermitFunction::HermitFunction(Eigen::Vector2d b, Eigen::Vector2d e, double angleb, double anglee)
{
    Eigen::MatrixXd P(2, 2);
    P << b.x(), e.x(), b.y(), e.y();
    m_P = P;

    angleb = angleb / 180.0f * 3.1415926f;
    anglee = anglee / 180.0f * 3.1415926f;
    Eigen::Vector2d p0 = Eigen::Vector2d(std::cosf(angleb), std::sinf(angleb));

    Eigen::Vector2d p1 = Eigen::Vector2d(std::cosf(anglee), std::sinf(anglee));
    Eigen::MatrixXd T(2, 2);
    T << p0.x(), p1.x(), p0.y(), p1.y();
    m_T = T;
}

void HermitFunction::HermitFunction2(Eigen::Vector2d b, Eigen::Vector2d e, double angleb, double anglee)
{
    Eigen::MatrixXd P(2, 2);
    P << b.x(), e.x(), b.y(), e.y();
    m_P = P;

    angleb = angleb / 180.0f * 3.1415926f;
    anglee = anglee / 180.0f * 3.1415926f;
    Eigen::Vector2d p0 = Eigen::Vector2d(std::cosf(angleb), std::sinf(angleb));

    Eigen::Vector2d p1 = Eigen::Vector2d(std::cosf(anglee), std::sinf(anglee));
    Eigen::MatrixXd T(2, 2);
    T << p0.x(), p1.x(), p0.y(), p1.y();
    m_T = T;
}

void HermitFunction::UsingpointWithtangent(Eigen::Vector2d b,
    Eigen::Vector2d p0,
    Eigen::Vector2d e,
    Eigen::Vector2d p1)
{
    Eigen::MatrixXd P(2, 2);
    P << b.x(), e.x(), b.y(), e.y();
    m_P = P;

    Eigen::MatrixXd T(2, 2);
    T << p0.x(), p1.x(), p0.y(), p1.y();
    m_T = T;
}

Eigen::Vector2d HermitFunction::getpoint(double t)
{
    int n = m_P.cols();
    int i = int(t * (n - 1));
    t = t * (n - 1) - i;

    Eigen::Vector2d point_on_curve = (2 * t * t * t - 3 * t * t + 1) * m_P.col(i)
        + (t * t * t - 2 * t * t + t) * m_T.col(i)
        + (-2 * t * t * t + 3 * t * t) * m_P.col(i + 1)
        + (t * t * t - t * t) * m_T.col(i + 1);

    return point_on_curve;
}

Eigen::Vector2d HermitFunction::getderivation(double t) 
{
    int n = m_P.cols();
    int i = int(t * (n - 1));
    t = t * (n - 1) - i;

    Eigen::Vector2d pointder_on_curve = (6 * t * t - 6 * t) * m_P.col(i)
        + (3 * t * t - 4 * t + 1) * m_T.col(i)
        + (-6 * t * t + 6 * t) * m_P.col(i + 1)
        + (3 * t * t - 2 * t) * m_T.col(i + 1);

    Eigen::Vector2d dP = 3 * (1 - t) * (1 - t) * (m_P.col(i + 1) - m_P.col(i))
        + 6 * (1 - t) * t * (m_T.col(i + 1) - m_T.col(i))
        + 3 * t * t * (m_P.col(i + 1) - m_P.col(i));

    return dP;
}

Eigen::Vector2d HermitFunction::GetPoint(double t) 
{
    Eigen::Vector2d P0_ = m_P.col(0);
    Eigen::Vector2d P1_ = m_P.col(1);

    Eigen::Vector2d T0_ = m_T.col(0);
    Eigen::Vector2d T1_ = m_T.col(1);
    double t2 = t * t;
    double t3 = t2 * t;
    double h1 = 2 * t3 - 3 * t2 + 1;
    double h2 = -2 * t3 + 3 * t2;
    double h3 = t3 - 2 * t2 + t;
    double h4 = t3 - t2;
    return h1 * P0_ + h2 * P1_ + h3 * T0_ + h4 * T1_;
}

std::vector<Eigen::Vector2d> HermitFunction::getPointsvec(int Numpoints) 
{
    double dert = 1.0f / (double)(Numpoints - 1);
    std::vector<Eigen::Vector2d> vec;
    for (double t = 0; t <= 1; t += dert) {
        Eigen::Vector2d point_on_curve = HermitFunction::getpoint(t);
        vec.push_back(point_on_curve);
    }
    return vec;
}

std::vector<Eigen::Vector2d> HermitFunction::getPointsvec(std::vector<Eigen::Vector2d> points,
    std::vector<Eigen::Vector2d> tangents,
    int Numpoints)
{
    double dert = 1.0f / (double)(Numpoints);
    std::vector<Eigen::Vector2d> vec;
    for (int i = 0; i < points.size() - 1; i++) {
        if (i % 2 == 0)
            UsingpointWithtangent(points[i], tangents[i], points[i + 1], tangents[i + 1]);
        else
            UsingpointWithtangent(points[i], tangents[i], points[i + 1], tangents[i + 1]);

        for (double t = 0; t <= 1; t += dert) {
            Eigen::Vector2d point_on_curve = HermitFunction::getpoint(t);
            vec.push_back(point_on_curve);
        }
    }

    return vec;
}

double HermitFunction::distance(const Eigen::Vector2d& point) 
{
    Eigen::Vector2d P0_ = m_P.col(0);
    Eigen::Vector2d P1_ = m_P.col(1);

    Eigen::Vector2d T0_ = m_T.col(0);
    Eigen::Vector2d T1_ = m_T.col(1);

    double epsilon = 1e-4;
    double t0 = 0;
    double t1 = 1;
    double t = 0.5;
    Eigen::Vector2d P = GetPoint(t);
    Eigen::Vector2d dP = 3 * (1 - t) * (1 - t) * (P1_ - P0_) + 6 * (1 - t) * t * (T1_ - T0_)
        + 3 * t * t * (P1_ - P0_);
    Eigen::Vector2d v = point - P;
    double dist = v.norm();
    while (dist > epsilon) {
        double a = dP.dot(dP);
        double b = 2 * v.dot(dP);
        double c = v.dot(v);
        double D = b * b - 4 * a * c;
        double t2 = (-b + sqrt(D)) / (2 * a);
        double t3 = (-b - sqrt(D)) / (2 * a);
        if (t2 >= 0 && t2 <= 1) {
            t = t2;
        }
        else if (t3 >= 0 && t3 <= 1) {
            t = t3;
        }
        else {
            t = (t2 + t3) / 2;
        }
        P = GetPoint(t);
        dP = 3 * (1 - t) * (1 - t) * (P1_ - P0_) + 6 * (1 - t) * t * (T1_ - T0_)
            + 3 * t * t * (P1_ - P0_);
        v = point - P;
        dist = v.norm();
    }
    return dist;
}


