#include "tool.h"

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
            nnPArray(i, j) = i * nelx + j + 1;
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
            edofMatNew(i, j) = dofVector(edofMat(i, j) - 1);
    edofMat = edofMatNew;
    int ndof = 2 * nnP; // Number of dofs

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

    //MatrixXf M1(3, 3);    // Column-major storage
    //M1 << 1, 2, 3,
    //	4, 5, 6,
    //	7, 8, 9;

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

    // Solve (remember to constrain one node)
    //Eigen::SparseMatrix<double> chi(ndof, 3);
    //Eigen::SparseMatrix<double> K_block = K.block(2, 2, ndof - 2, ndof - 2);
    //chi.block(2, 0, ndof - 2, 3) = K_block.llt().solve(F.block(2, 0, ndof - 2, 3));

    SparseLU<SparseMatrix<double>> solver;
    solver.compute(K.transpose().block(2,
        2,
        ndof - 2,
        ndof - 2)); // 求解K中3到ndof的部分，注意要去除边界上的自由度

    SparseMatrix<double> chi;

    chi = solver.solve(F.block(2, 0, ndof - 2, 3)); // 求解线性方程组
    //std::cout << chi << '\n';

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

    //std::cout << "ke3_5end_3_5end:\n" << ke3_5end_3_5end << "++++++++++++++++++++++\n";
    //std::cout << "fe3_5end\n" << fe3_5end << "++++++++++++++++++++++\n";

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

            //std::cout << "sumLambda:\n" << sumLambda << "++++++++++++++++++++++\n";
            //std::cout << "LambdaMatrix:\n" << LambdaMatrix << "++++++++++++++++++++++\n";

            //std::cout << "LambdaVec.*sumLambda:\n" << LambdaMatrix.cwiseProduct(sumLambda) << "+++++++++++++++++++\n";
            //std::cout << "MuVec.*sumMu:\n" << MuMatrix.cwiseProduct(sumMu) << "+++++++++++++++++++\n";

            CH(i, j) = 1.0 / cellVolume
                * (LambdaMatrix.cwiseProduct(sumLambda) + MuMatrix.cwiseProduct(sumMu)).sum();
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
            nnPArray(i, j) = i * nelx + j + 1;
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
            edofMatNew(i, j) = dofVector(edofMat(i, j) - 1);
    edofMat = edofMatNew;
    int ndof = 2 * nnP; // Number of dofs

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

    //MatrixXf M1(3, 3);    // Column-major storage
    //M1 << 1, 2, 3,
    //	4, 5, 6,
    //	7, 8, 9;

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

    // Solve (remember to constrain one node)
    //Eigen::SparseMatrix<double> chi(ndof, 3);
    //Eigen::SparseMatrix<double> K_block = K.block(2, 2, ndof - 2, ndof - 2);
    //chi.block(2, 0, ndof - 2, 3) = K_block.llt().solve(F.block(2, 0, ndof - 2, 3));

    SparseLU<SparseMatrix<double>> solver;
    solver.compute(K.transpose().block(2,
        2,
        ndof - 2,
        ndof - 2)); // 求解K中3到ndof的部分，注意要去除边界上的自由度

    SparseMatrix<double> chi;

    chi = solver.solve(F.block(2, 0, ndof - 2, 3)); // 求解线性方程组
    //std::cout << chi << '\n';

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

    //std::cout << "ke3_5end_3_5end:\n" << ke3_5end_3_5end << "++++++++++++++++++++++\n";
    //std::cout << "fe3_5end\n" << fe3_5end << "++++++++++++++++++++++\n";

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

            //std::cout << "sumLambda:\n" << sumLambda << "++++++++++++++++++++++\n";
            //std::cout << "LambdaMatrix:\n" << LambdaMatrix << "++++++++++++++++++++++\n";

            //std::cout << "LambdaVec.*sumLambda:\n" << LambdaMatrix.cwiseProduct(sumLambda) << "+++++++++++++++++++\n";
            //std::cout << "MuVec.*sumMu:\n" << MuMatrix.cwiseProduct(sumMu) << "+++++++++++++++++++\n";

            CH(i, j) = 1.0 / cellVolume
                * (LambdaMatrix.cwiseProduct(sumLambda) + MuMatrix.cwiseProduct(sumMu)).sum();
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