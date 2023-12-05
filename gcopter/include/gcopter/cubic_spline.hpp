/*
 * @Author: Shirley
 * @Date: 2022-09-20 12:32:27
 * @LastEditors: Shirley
 * @LastEditTime: 2022-10-08 22:11:13
 * @Description: file content
 */
#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "cubic_curve.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <vector>
// 三次样条曲线
// 三次样条曲线
// 三次样条曲线

namespace cubic_spline
{
    // N*N注意注意注意：这里的N =  传入的参数N-1！！！！！！！！！！！！！！！！！！！
    // The banded system class is used for solving
    // banded linear system Ax=b efficiently.
    // A is an N*N band matrix with lower band width lowerBw
    // and upper band width upperBw.
    // Banded LU factorization has O(N) time complexity.
    // LU分解
    class BandedSystem
    {
    public:
        // The size of A, as well as the lower/upper
        // banded width p/q are needed
        inline void create(const int &n, const int &p, const int &q)
        {
            // In case of re-creating before destroying
            destroy();
            N = n;
            lowerBw = p;
            upperBw = q;
            int actualSize = N * (lowerBw + upperBw + 1);
            ptrData = new double[actualSize];
            std::fill_n(ptrData, actualSize, 0.0);
            return;
        }
        inline void destroy()
        {
            // 经典的用法delete之后再让其指向nullptr
            if (ptrData != nullptr)
            {
                delete[] ptrData;
                ptrData = nullptr;
            }
            return;
        }

    private:
        int N;
        int lowerBw;
        int upperBw;
        // 将矩阵存放在一个一维数组中 size:N*N!
        // Compulsory nullptr initialization here
        double *ptrData = nullptr;

    public:
        // Reset the matrix to zero
        inline void reset(void)
        {
            std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
            return;
        }
        // 该类实现了()运算符重载，用于获取矩阵中(i,j)位置处的元素
        // The band matrix is stored as suggested in "Matrix Computation"
        inline const double &operator()(const int &i, const int &j) const
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        inline double &operator()(const int &i, const int &j)
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        // This function conducts banded LU factorization in place
        // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
        inline void factorizeLU()
        {
            int iM, jM;
            double cVl;
            for (int k = 0; k <= N - 2; ++k)
            {
                iM = std::min(k + lowerBw, N - 1);
                cVl = operator()(k, k);
                for (int i = k + 1; i <= iM; ++i)
                {
                    if (operator()(i, k) != 0.0)
                    {
                        operator()(i, k) /= cVl;
                    }
                }
                jM = std::min(k + upperBw, N - 1);
                for (int j = k + 1; j <= jM; ++j)
                {
                    cVl = operator()(k, j);
                    if (cVl != 0.0)
                    {
                        for (int i = k + 1; i <= iM; ++i)
                        {
                            if (operator()(i, k) != 0.0)
                            {
                                operator()(i, j) -= operator()(i, k) * cVl;
                            }
                        }
                    }
                }
            }
            return;
        }

        // This function solves Ax=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solve(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; ++j)
            {
                iM = std::min(j + lowerBw, N - 1);
                for (int i = j + 1; i <= iM; ++i)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j)
            {
                b.row(j) /= operator()(j, j);
                iM = std::max(0, j - upperBw);
                for (int i = iM; i <= j - 1; ++i)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            return;
        }

        // This function solves ATx=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solveAdj(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; ++j)
            {
                b.row(j) /= operator()(j, j);
                iM = std::min(j + upperBw, N - 1);
                for (int i = j + 1; i <= iM; ++i)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j)
            {
                iM = std::max(0, j - lowerBw);
                for (int i = iM; i <= j - 1; ++i)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
        }
    };
    // 三次样条曲线
    class CubicSpline
    {
    public:
        CubicSpline() = default;
        ~CubicSpline() { A.destroy(); }

    private:
        int N;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        BandedSystem A; // 里面的矩阵是n-1阶方阵，但是是用一维数组表示和存储的！
        Eigen::MatrixX2d b; // n-1行
        Eigen::MatrixX2d allP;

        Eigen::MatrixXd dB_dx; // 干啥的？？
        Eigen::MatrixXd dxx_dx;
        std::vector<Eigen::Matrix<double, 2, 4>> cMats; // 系数矩阵，一个维度一行

    public:
        inline Eigen::Vector2d getHeadP() const { return headP; }
        inline Eigen::Vector2d getTailP() const { return tailP; }
        // 设置端点
        inline void
        setConditions(const Eigen::Vector2d &headPos,
                      const Eigen::Vector2d &tailPos,
                      const int &pieceNum)
        {
            // TODO
            headP = headPos;
            tailP = tailPos;
            N = pieceNum;
            // A：N-1方阵 b:N-1行
            A.create(N - 1, 2, 2);

            b.resize(N - 1, 2);

            allP.resize(N + 1, 2);
            allP.row(0) = headP;
            allP.row(N) = tailP;

            dB_dx.resize(N - 1, N - 1);
            dB_dx.setZero();
            for (int i = 0; i < N - 2; ++i)
            {
                dB_dx(i, i + 1) = 3;
                dB_dx(i + 1, i) = -3;
            }

            dxx_dx.resize(N, N - 1);
            dxx_dx.setZero();
            for (int i = 0; i < N - 1; ++i)
            {
                dxx_dx(i, i) = -1;
                dxx_dx(i + 1, i) = 1;
            }
            cMats.resize(N);

            return;
        }
        // 设置内点
        inline void setInnerPoints(const Eigen::Ref<const Eigen::Matrix2Xd> &inPs)
        {
            //  TODO
            A.reset();
            b.setZero();
            for (int i = 0; i < N - 1; i++)
            {
                allP.row(i + 1) = inPs.col(i).transpose();
            }
            for (int i = 0; i < N - 2; i++)
            {
                b.row(i) = 3 * (allP.row(i + 2) - allP.row(i));

                A(i, i) = 4;
                A(i + 1, i) = 1;
                A(i, i + 1) = 1;
            }

            b.row(N - 2) = 3 * (allP.row(N) - allP.row(N - 2));
            A(N - 2, N - 2) = 4;

            A.factorizeLU();
            A.solve(b);
            A.solve(dB_dx); // dD_dx is stored in dB_dx, (n-1)*(n-1)

            Eigen::Matrix<double, 2, 4> cMat;
            cMat.col(3) = allP.row(0).transpose();
            cMat.col(2).setZero();
            cMat.col(1) = (3 * (allP.row(1) - allP.row(0)) - b.row(0)).transpose();
            cMat.col(0) = (2 * (allP.row(0) - allP.row(1)) + b.row(0)).transpose();
            cMats[0] = cMat;

            for (int i = 1; i < N - 1; ++i)
            {
                cMat.col(3) = allP.row(i).transpose();
                cMat.col(2) = b.row(i - 1).transpose();
                cMat.col(1) = (3 * (allP.row(i + 1) - allP.row(i)) - b.row(i) - 2 * b.row(i - 1)).transpose();
                cMat.col(0) = (2 * (allP.row(i) - allP.row(i + 1)) + b.row(i - 1) + b.row(i)).transpose();
                cMats[i] = cMat;
            }
            cMat.col(3) = allP.row(N - 1).transpose();
            cMat.col(2).setZero();
            cMat.col(1) = (3 * (allP.row(N) - allP.row(N - 1)) - 2 * b.row(N - 2)).transpose();
            cMat.col(0) = (2 * (allP.row(N - 1) - allP.row(N)) + b.row(N - 2)).transpose();
            cMats[N - 1] = cMat;
            return;
        }
        // 设置内点，2d中间点使用列向量表示
        inline void setInnerPointsV(const Eigen::Ref<const Eigen::VectorXd> &inPs)
        {
            // TODO
            A.reset();
            b.setZero();
            for (int i = 0; i < N - 1; i++)
            {
                // 上下两个函数，只有这个地方不一样，其他地方均相同。
                // 2d,依此为x、y
                allP(i + 1, 0) = inPs(i);
                allP(i + 1, 1) = inPs(i + N - 1);
            }
            for (int i = 0; i < N - 2; i++)
            {
                b.row(i) = 3 * (allP.row(i + 2) - allP.row(i));

                A(i, i) = 4;
                A(i + 1, i) = 1;
                A(i, i + 1) = 1;
            }
            b.row(N - 2) = 3 * (allP.row(N) - allP.row(N - 2));
            A(N - 2, N - 2) = 4;
            // 三角分解LU分解，详细见东北大学矩阵分析教材
            A.factorizeLU();
            A.solve(b);
            A.solve(dB_dx); // dD_dx is stored in dB_dx, (n-1)*(n-1)

            // 求解完线性方程组之后，进行三次样条曲线的参数存储
            Eigen::Matrix<double, 2, 4> cMat;
            cMat.col(3) = allP.row(0).transpose();
            cMat.col(2).setZero();
            cMat.col(1) = (3 * (allP.row(1) - allP.row(0)) - b.row(0)).transpose();
            cMat.col(0) = (2 * (allP.row(0) - allP.row(1)) + b.row(0)).transpose();
            cMats[0] = cMat;
            for (int i = 1; i < N - 1; ++i)
            {
                cMat.col(3) = allP.row(i).transpose();
                cMat.col(2) = b.row(i - 1).transpose();
                cMat.col(1) = (3 * (allP.row(i + 1) - allP.row(i)) - b.row(i) - 2 * b.row(i - 1)).transpose();
                cMat.col(0) = (2 * (allP.row(i) - allP.row(i + 1)) + b.row(i - 1) + b.row(i)).transpose();

                cMats[i] = cMat;
            }
            cMat.col(3) = allP.row(N - 1).transpose();
            cMat.col(2).setZero();
            cMat.col(1) = (3 * (allP.row(N) - allP.row(N - 1)) - 2 * b.row(N - 2)).transpose();
            cMat.col(0) = (2 * (allP.row(N - 1) - allP.row(N)) + b.row(N - 2)).transpose();
            cMats[N - 1] = cMat;
            return;
        }
        inline void getCurve(CubicCurve &curve) const
        {
            // TODO
            curve.clear();
            curve.reserve(N);
            std::vector<double> durs(N, 1.0);
            auto c = CubicCurve(durs, cMats);
            curve = std::move(c);
            return;
        }
        // 得到能量函数，牛
        inline void getStretchEnergy(double &energy) const
        {
            // TODO
            energy = 0;
            for (int i = 0; i < N; ++i)
            {
                energy += 4 * cMats[i].col(1).dot(cMats[i].col(1)) +
                          12 * cMats[i].col(1).dot(cMats[i].col(0)) +
                          12 * cMats[i].col(0).dot(cMats[i].col(0));
            }
            return;
        }

        inline const Eigen::MatrixX2d &getCoeffs(void) const
        {
            return b;
        }
        // 设置梯度
        inline void getGrad(Eigen::Ref<Eigen::Matrix2Xd> gradByPoints) const
        {
            // TODO
            gradByPoints.setZero();
            Eigen::MatrixXd ddi_dx = (2 * dxx_dx.row(0) + dB_dx.row(0)).transpose();
            Eigen::MatrixXd dci_dx = (-3 * dxx_dx.row(0) - dB_dx.row(0)).transpose();
            gradByPoints += (24 * cMats[0].col(0) * ddi_dx.transpose() + 12 * cMats[0].col(0) * dci_dx.transpose() +
                             12 * cMats[0].col(1) * ddi_dx.transpose() + 8 * cMats[0].col(1) * dci_dx.transpose());
            for (int i = 1; i < N - 1; ++i)
            {
                ddi_dx = (2 * dxx_dx.row(i) + dB_dx.row(i - 1) + dB_dx.row(i)).transpose();
                dci_dx = (-3 * dxx_dx.row(i) - 2 * dB_dx.row(i - 1) - dB_dx.row(i)).transpose();
                gradByPoints += (24 * cMats[i].col(0) * ddi_dx.transpose() + 12 * cMats[i].col(0) * dci_dx.transpose() +
                                 12 * cMats[i].col(1) * ddi_dx.transpose() + 8 * cMats[i].col(1) * dci_dx.transpose());
            }

            ddi_dx = (2 * dxx_dx.row(N - 1) + dB_dx.row(N - 2)).transpose();
            dci_dx = (-3 * dxx_dx.row(N - 1) - 2 * dB_dx.row(N - 2)).transpose();
            gradByPoints += (24 * cMats[N - 1].col(0) * ddi_dx.transpose() + 12 * cMats[N - 1].col(0) * dci_dx.transpose() +
                             12 * cMats[N - 1].col(1) * ddi_dx.transpose() + 8 * cMats[N - 1].col(1) * dci_dx.transpose());
        }
    };
}

#endif
