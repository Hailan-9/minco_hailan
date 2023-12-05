/*
 * @Author: Wang Zhepei
 * @Date: 2022-08-05 17:46:25
 * @LastEditors: “Shirley” “todo@163.com”
 * @LastEditTime: 2022-10-09 10:22:39
 * @Description: file content
 */
#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include "cubic_spline.hpp"
#include "lbfgs.hpp"
#include "sdqp.hpp"
#include <Eigen/Eigen>

#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>

namespace path_smoother
{
    class PathSmoother
    {
    private:
        cubic_spline::CubicSpline cubSpline;

        int pieceN;
        Eigen::Matrix3Xd diskObstacles;
        std::vector<Eigen::MatrixX3d> hPolysObstacles;
        double penaltyWeight;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        Eigen::Matrix2Xd points;
        Eigen::Matrix2Xd gradByPoints;

        lbfgs::lbfgs_parameter_t lbfgs_params;

    private:
        static inline double costFunction(void *ptr,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &g)
        {
            // TODO
            PathSmoother &obj = *(PathSmoother *)ptr;
            g.setZero();
            double energy_cost, potential_cost = 0;
            double x_num = obj.pieceN - 1;
            obj.cubSpline.setInnerPointsV(x);
            obj.cubSpline.getStretchEnergy(energy_cost);
            obj.cubSpline.getGrad(obj.gradByPoints);
            Eigen::Map<const Eigen::MatrixXd> xx(x.data(), x_num, 2);
            Eigen::Map<Eigen::MatrixX2d> gg(g.data(), x_num, 2);

            Eigen::Matrix<double, 2, 2> Q = 2 * Eigen::Matrix2d::Identity();
            double poly_num = obj.hPolysObstacles.size();
            for (int i = 0; i < x_num; ++i)
            {
                /** robot position pr = (xr, yr) ^T,  the nearest position in obstacle is p = (x, y)^T
                 * min(0.5 * p' 2*Indentity p + ( -2) pr' p )
                 * s.t. poly[ : , 0:1] p <= poly[ : , 2]
                 *
                 * sdqp:
                 *  minimize     0.5 x' Q x + c' x
                 * subject to       A x <= b
                 * Q must be positive definite
                 * */
                Eigen::Matrix<double, 2, 1> x_obs;
                x_obs.setZero();
                Eigen::Matrix<double, 2, 1> x_cur = xx.row(i).transpose();
                Eigen::Matrix<double, 2, 1> c = -2 * x_cur;
                for (int j = 0; j < poly_num; ++j)
                {
                    Eigen::Matrix<double, -1, 2> A = obj.hPolysObstacles.at(j).leftCols(2);
                    Eigen::Matrix<double, -1, 1> b = -obj.hPolysObstacles.at(j).rightCols(1);
                    Eigen::Matrix<double, -1, 1> coeff_to_planes = A * x_cur - b;
                    if (coeff_to_planes.maxCoeff() > 0) // out obstacle
                    {
                        sdqp::sdqp<2>(Q, c, A, b, x_obs);
                        Eigen::Matrix<double, 2, 1> dist_vec = x_cur - x_obs;
                        double dist = dist_vec.norm();
                        if (dist < 1.0)
                        {
                            double dist_cost = 1000 * (1 - dist);
                            potential_cost += dist_cost;
                            Eigen::Matrix<double, 1, 2> cur_g = -1000 / dist * dist_vec.transpose();
                            gg.row(i) += cur_g;
                        }
                    }
                    else // in obstacle
                    {
                        /***method 1***/
                        Eigen::Matrix<double, -1, 1> planes_vector_normals = A.rowwise().norm(); // sqrt(a^2 + b^2)
                        Eigen::Matrix<double, -1, 1> tmp = coeff_to_planes.array().abs();
                        Eigen::Matrix<double, -1, 1> dist_to_planes = tmp.cwiseQuotient(planes_vector_normals);
                        int r, co;
                        double dist_cost = 1000 * (dist_to_planes.minCoeff(&r, &co) + 1);
                        potential_cost += dist_cost;
                        Eigen::Matrix<double, 1, 2> cur_g = -1000 * A.row(r) / planes_vector_normals(r);
                        gg.row(i) += cur_g;

                        /***method 2***/
                        // double a = A(r, 0), bb = A(r, 1), c = -b(r, 0);
                        // double x_nearest = (bb * (bb * x_cur(0) - a * x_cur(1)) - a * c) / (a * a + bb * bb);
                        // double y_nearest = (a * (-bb * x_cur(0) + a * x_cur(1)) - bb * c) / (a * a + bb * bb);
                        // Eigen::Matrix<double, 2, 1> nearest_vec(x_cur(0) - x_nearest, x_cur(1) - y_nearest);
                        // double nearest_dist = nearest_vec.norm();
                        // Eigen::Matrix<double, 1, 2> nearest_gg_cur = nearest_vec.transpose() / nearest_dist;

                        // potential_cost += 1000 * (nearest_dist + 1);
                        // gg.row(i) += 1000 * nearest_gg_cur;
                    }
                }
            }
            double cost = obj.penaltyWeight * energy_cost + potential_cost;

            gg = gg + obj.penaltyWeight * obj.gradByPoints.transpose();

            return cost;
        }

    public:
        /*Curcle obstacles*/
        inline bool setup(const Eigen::Vector2d &initialP,
                          const Eigen::Vector2d &terminalP,
                          const int &pieceNum,
                          const Eigen::Matrix3Xd &diskObs,
                          const double penaWeight)
        {
            pieceN = pieceNum;
            diskObstacles = diskObs;
            penaltyWeight = penaWeight;
            headP = initialP;
            tailP = terminalP;

            cubSpline.setConditions(headP, tailP, pieceN);
            points.resize(2, pieceN - 1);
            gradByPoints.resize(2, pieceN - 1);
            return true;
        }
        /* Polytope obstacles 多胞形*/
        inline bool setup(const Eigen::Vector2d &initialP,
                          const Eigen::Vector2d &terminalP,
                          const int &pieceNum,
                          std::vector<Eigen::MatrixX3d> &hPolys,
                          const double penaWeight)
        {
            std::cout << "setup" << std::endl;
            pieceN = pieceNum;

            hPolysObstacles = hPolys;
            penaltyWeight = penaWeight;
            headP = initialP;
            tailP = terminalP;

            cubSpline.setConditions(headP, tailP, pieceN);
            points.resize(2, pieceN - 1);
            gradByPoints.resize(2, pieceN - 1);
            return true;
        }

        inline double optimize(CubicCurve &curve,
                               const Eigen::Matrix2Xd &iniInPs,
                               const double &relCostTol)
        {
            // TODO
            cubSpline.setInnerPoints(iniInPs);

            int k = pieceN - 1;
            // 存储内点的数组，2d的航迹中间点使用一个列向量表示，这些即决策变量！！！
            Eigen::VectorXd x(k * 2);

            for (int i = 0; i < k; ++i)
            {
                x(i) = iniInPs(0, i);
                x(i + k) = iniInPs(1, i);
            }

            double f;
            lbfgs::lbfgs_parameter_t lbfgs_params;
            lbfgs_params.mem_size = 256;
            lbfgs_params.past = 3;
            lbfgs_params.min_step = 1.0e-32;
            lbfgs_params.max_step = 100;
            lbfgs_params.g_epsilon = 1e-6;
            lbfgs_params.delta = relCostTol;
            lbfgs_params.max_iterations = 10000;
            lbfgs_params.use_alm = false;

            int ret = lbfgs::lbfgs_optimize(x,
                                            f,
                                            &PathSmoother::costFunction,
                                            nullptr,
                                            nullptr,
                                            this,
                                            lbfgs_params);
            if (ret < 0)
            {
                std::cout << "\033[32m" << lbfgs::lbfgs_strerror(ret) << "\033[0m" << std::endl;
                return ret;
            }
            double minCost;
            if (ret >= 0)
            {
                minCost = f;
                cubSpline.setInnerPointsV(x);
                cubSpline.getCurve(curve);
                std::cout << "[path_smoother] success minCost = " << minCost << std::endl;
            }

            return minCost;
        }
    };
}

#endif
