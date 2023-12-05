/*
 * @Author: Shirley
 * @Date: 2022-10-06 15:07:08
 * @LastEditors: “Shirley” “todo@163.com”
 * @LastEditTime: 2022-10-09 18:33:32
 * @Description: file content
 */

#ifndef TOPP_HPP
#define PTOPP_HPP

#include "lbfgs.hpp"
#include <Eigen/Eigen>

#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>
namespace topp
{
    class Topp
    {
    public:
        Topp() = default;

        Topp(const CubicCurve &curve, const double &resolution)
        {
            // N段轨迹
            const int N = curve.getPieceNum();
            K = N * resolution;
            total_t = curve.getTotalDuration();
            dt = total_t / K;
            // 锥标准
            cone_criterion.resize(3, K * 2 + 1);
            projection_on_symm_cone.resize(3, K * 2 + 1);
            mu.resize(3, K * 2 + 1);
            mu.setZero();
            mu_rho.resize(3, K * 2 + 1);

            lambda_rho.resize(K + 2);
            lambda.resize(K + 2);
            lambda.setZero();
            equal_criterion.resize(K + 2);
            sigma_rho.resize(6 * K + 2);
            sigma.resize(6 * K + 2);
            sigma.setZero();
            inequal_criterion.resize(6 * K + 2);

            /************ ******************************
             * calculate s, q, q' , q''
            s is not uniform distribution
            s = last_s + (v_last + v_cur)/2 * dt
            v = a1 + 2*a1*t + 3*a3*t
            q = p
            q' = p' * dt/ds
            qx'' = (px'' *  py' - py'' * px')  * py' * (dt/ds)^4
            *******************************************/
            s.resize(K + 1);
            v_list.resize(K + 1);
            q.resize(2, K + 1);
            dq.resize(2, K + 1);
            ddq.resize(2, K);
            dq_norm2.resize(K + 1);
            ax.resize(K);
            ay.resize(K);

            Eigen::Vector2d q_vec, v_vec, a_vec;

            for (int k = 0; k < K + 1; ++k)
            {
                double t = dt * k;
                curve.getPosVelAcc(t, q_vec, v_vec, a_vec);
                int count = 1000;
                double v_norm = v_vec.norm();
                while (v_norm < 1e-10 && --count)
                {
                    if (t < total_t / 2)
                        t += 1e-3;
                    else
                        t -= 1e-3;
                    curve.getPosVelAcc(t, q_vec, v_vec, a_vec);
                    v_norm = v_vec.norm();
                }
                v_list(k) = v_norm;
                double dt_ds = 1 / v_list(k);
                dq.col(k) = v_vec * dt_ds;
                // q.col(k) = q_vec;
                // double dt_ds4 = dt_ds * dt_ds * dt_ds * dt_ds;
                // ddq(0, k) = (a_vec(0) * v_vec(1) - a_vec(1) * v_vec(0)) * v_vec(1) * dt_ds4;
                // ddq(1, k) = (a_vec(1) * v_vec(0) - a_vec(0) * v_vec(1)) * v_vec(0) * dt_ds4;
            }
            dq_norm2 = dq.colwise().squaredNorm(); // equals 1, in fact it is useless
            s(0) = 0;
            for (int k = 1; k < K + 1; ++k)
            {
                s(k) = s(k - 1) + (v_list(k) + v_list(k - 1)) / 2 * dt;
            }
            for (int k = 0; k < K; ++k)
            {
                ddq.col(k) = (dq.col(k + 1) - dq.col(k)) / (s(k + 1) - s(k) + 1e-10);
            }
        }
        Eigen::VectorXd ax, ay, s;

    private:
        lbfgs::lbfgs_parameter_t lbfgs_params;
        double dt, total_t;
        int K;
        double v_max_, a_max_, v_start_, v_end_;

        double rho;
        const double gamma_ = 1., beta_ = 1e3;
        Eigen::Matrix3Xd mu_rho, mu, cone_criterion, projection_on_symm_cone;
        Eigen::Matrix2Xd q, dq, ddq;
        Eigen::VectorXd dq_norm2;
        Eigen::VectorXd lambda_rho, lambda, equal_criterion, sigma_rho, sigma, inequal_criterion;
        Eigen::VectorXd v_list; // ds/dt

    public:
        /*project on symm cone*/
        inline void printInfo()
        {
            std::cout << "dq = \n"
                      << dq << std::endl;
            std::cout << "ddq = \n"
                      << ddq << std::endl;
            std::cout << "dq_norm2 = \n"
                      << dq_norm2 << std::endl;
            std::cout << "s = \n"
                      << s << std::endl;
        }

        inline Eigen::Vector3d projectVectorOnCone(Eigen::Vector3d u)
        {
            Eigen::Vector3d projection;
            int dim = u.size() - 1;
            double u0 = u(0);
            Eigen::Map<Eigen::VectorXd> u1(u.data() + 1, dim);
            double u1_norm = u1.norm();
            if (u0 <= -u1_norm)
            {
                projection.setZero();
            }
            else if (u0 >= u1_norm)
            {
                projection = u;
            }
            else
            {
                projection(0) = (u0 + u1_norm) / 2;
                projection.segment(1, dim) = u1 * (u0 + u1_norm) / (2 * u1_norm);
            }
            return projection;
        }

        /* cost function for task2*/
        static inline double costFunction(void *ptr,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &g)
        {
            Topp &obj = *(Topp *)ptr;
            g.setZero();
            double cost = 0;
            const int K = obj.K;
            const double rho = obj.rho;
            Eigen::Map<const Eigen::VectorXd> bk(x.data(), K + 1);
            Eigen::Map<const Eigen::VectorXd> ak(x.data() + K + 1, K);
            Eigen::Map<const Eigen::VectorXd> ck(x.data() + 2 * K + 1, K + 1);
            Eigen::Map<const Eigen::VectorXd> dk(x.data() + 3 * K + 2, K);
            Eigen::Map<Eigen::VectorXd> gbk(g.data(), K + 1);
            Eigen::Map<Eigen::VectorXd> gak(g.data() + K + 1, K);
            Eigen::Map<Eigen::VectorXd> gck(g.data() + 2 * K + 1, K + 1);
            Eigen::Map<Eigen::VectorXd> gdk(g.data() + 3 * K + 2, K);
            /*** conic alm***/
            Eigen::Vector3d u, proj_u;
            gck(0) = 0.0;
            for (int k = 0; k < K; ++k)
            {
                u(0) = ck(k + 1) + ck(k) + dk(k);
                u(1) = 2;
                u(2) = ck(k + 1) + ck(k) - dk(k);
                u = obj.mu_rho.col(k) - u;
                proj_u = obj.projectVectorOnCone(u);

                gck(k) += -rho * (proj_u(0) + proj_u(2));
                gck(k + 1) += -rho * (proj_u(0) + proj_u(2));
                gdk(k) += -rho * (proj_u(0) - proj_u(2));
                cost += proj_u.dot(proj_u);
                obj.projection_on_symm_cone.col(k) = proj_u;
            }

            for (int k = 0; k < K + 1; ++k)
            {
                u(0) = bk(k) + 1;
                u(1) = 2 * ck(k);
                u(2) = bk(k) - 1;
                u = obj.mu_rho.col(K + k) - u;
                proj_u = obj.projectVectorOnCone(u);
                obj.projection_on_symm_cone.col(K + k) = proj_u;
                gbk(k) += -rho * (proj_u(0) + proj_u(2));
                gck(k) += -rho * 2 * proj_u(1);
                cost += proj_u.dot(proj_u);
            }

            /*** equality constraint ***/
            {
                double t;
                Eigen::Vector3d u1, G(0., 1., -1.), gt;
                for (int k = 0; k < K; ++k)
                {
                    u1(0) = ak(k);
                    u1(1) = bk(k);
                    u1(2) = bk(k + 1);
                    G(0) = 2 * (obj.s(k + 1) - obj.s(k));
                    obj.equal_criterion(k) = u1.dot(G);
                    t = obj.equal_criterion(k) + obj.lambda_rho(k);
                    cost += t * t;
                    gt = rho * t * G;
                    gak(k) += gt(0);
                    gbk(k) += gt(1);
                    gbk(k + 1) += gt(2);
                }
                obj.equal_criterion(K) = obj.dq_norm2(0) * bk(0) - obj.v_start_ * obj.v_start_;
                t = obj.equal_criterion(K) + obj.lambda_rho(K);
                gbk(0) += rho * t * obj.dq_norm2(0);
                cost += t * t;
                obj.equal_criterion(K + 1) = obj.dq_norm2(K) * bk(K) - obj.v_end_ * obj.v_end_;
                t = obj.equal_criterion(K + 1) + obj.lambda_rho(K + 1);
                gbk(K) += rho * t * obj.dq_norm2(K);
                cost += t * t;
            }

            /*** inequality constraint ***/
            {
                obj.inequal_criterion = -obj.sigma_rho;
                obj.sigma.setZero();
                //-bk<=0
                double t;
                for (int k = 0; k < K + 1; ++k)
                {
                    t = -bk(k) + obj.sigma_rho(k);
                    if (t > 0)
                    {
                        gbk(k) += rho * t * (-1);
                        cost += t * t;
                        obj.sigma(k) = t * rho;
                        obj.inequal_criterion(k) = t - obj.sigma_rho(k);
                    }
                }

                // dq2 * bk <= v_max^2
                double dq2;
                for (int k = 0; k < K + 1; ++k)
                {
                    dq2 = obj.dq_norm2(k);
                    t = dq2 * bk(k) - obj.v_max_ * obj.v_max_ + obj.sigma_rho(K + 1 + k);
                    if (t > 0)
                    {
                        gbk(k) += rho * t * dq2;
                        cost += t * t;
                        obj.sigma(K + 1 + k) = t * rho;
                        obj.inequal_criterion(K + 1 + k) = t - obj.sigma_rho(K + 1 + k);
                    }
                }

                // v limit and a limit
                for (int k = 0; k < K; ++k)
                {
                    obj.ax(k) = bk(k) * obj.ddq(0, k) + obj.dq(0, k) * ak(k);
                    obj.ay(k) = bk(k) * obj.ddq(1, k) + obj.dq(1, k) * ak(k);
                }
                for (int k = 0; k < K; ++k)
                {
                    // x: ddq * bk + dq*ak<= a_max
                    t = obj.ax(k) - obj.a_max_ + obj.sigma_rho(2 * K + 2 + k);
                    if (t > 0)
                    {
                        gbk(k) += rho * t * obj.ddq(0, k);
                        gak(k) += rho * t * obj.dq(0, k);
                        cost += t * t;
                        obj.sigma(2 * K + 2 + k) = t * rho;
                        obj.inequal_criterion(2 * K + 2 + k) = obj.ax(k) - obj.a_max_;
                    }

                    // y: ddq * bk + dq*ak<= a_max
                    t = obj.ay(k) - obj.a_max_ + obj.sigma_rho(3 * K + 2 + k);
                    if (t > 0)
                    {
                        gbk(k) += rho * t * obj.ddq(1, k);
                        gak(k) += rho * t * obj.dq(1, k);
                        cost += t * t;
                        obj.sigma(3 * K + 2 + k) = t * rho;
                        obj.inequal_criterion(3 * K + 2 + k) = obj.ay(k) - obj.a_max_;
                    }

                    t = -obj.ax(k) - obj.a_max_ + obj.sigma_rho(4 * K + 2 + k);
                    if (t > 0)
                    {
                        gbk(k) += rho * t * (-obj.ddq(0, k));
                        gak(k) += rho * t * (-obj.dq(0, k));
                        cost += t * t;
                        obj.sigma(4 * K + 2 + k) = t * rho;
                        obj.inequal_criterion(4 * K + 2 + k) = -obj.ax(k) - obj.a_max_;
                    }

                    t = -obj.ay(k) - obj.a_max_ + obj.sigma_rho(5 * K + 2 + k);
                    if (t > 0)
                    {
                        gbk(k) += rho * t * (-obj.ddq(1, k));
                        gak(k) += rho * t * (-obj.dq(1, k));
                        cost += t * t;
                        obj.sigma(5 * K + 2 + k) = t * rho;
                        obj.inequal_criterion(5 * K + 2 + k) = -obj.ay(k) - obj.a_max_;
                    }
                }
            }
            /****optimal function****/
            double f = 2 * (obj.s.segment(1, K) - obj.s.segment(0, K)).dot(dk);
            gdk += 2 * (obj.s.segment(1, K) - obj.s.segment(0, K));

            return f + cost * rho / 2;
        }

        /* generate speed and acc for project task2*/
        inline double optimize(Eigen::VectorXd &bs,
                               Eigen::VectorXd &as,
                               const double &v_max,
                               const double &a_max,
                               const double &v_start,
                               const double &v_end,
                               const double &relCostTol)
        {
            std::cout << "topp optimize start" << std::endl;
            v_max_ = v_max;
            a_max_ = a_max;
            v_start_ = v_start;
            v_end_ = v_end;

            /*******************
             * bk : 0 ~ K (K+1)
             * ak : 0 ~ K-1 (K)
             * ck : 0 ~ K (K+1)
             * dk : 0 ~ K-1 (K)
             * ******************
             * x :
             * bk(0,K),
             * ak(K+1, 2*K),
             * ck(2*K+1,3*K+1),
             * dk(3*K+2,4*K+1),
             * *****************/
            Eigen::VectorXd x(4 * K + 2);
            // set init x
            for (int i = 0; i < K + 1; ++i)
            {
                x(i) = (v_end_ * v_end_ - v_start_ * v_start_) * i / K + v_start_ * v_start_; // bk
                x(2 * K + 1 + i) = sqrt(x(i));                                                // ck
            }
            double a_init = (v_end_ - v_start_) / total_t;
            for (int i = 0; i < K; ++i)
            {
                x(K + 1 + i) = a_init;                                              // ak
                x(3 * K + 2 + i) = 1.0 / (x(2 * K + 1 + i) + x(2 * K + 1 + i + 1)); // dk
            }

            /*********************ALM optimize***********************/
            lbfgs::lbfgs_parameter_t lbfgs_params;
            lbfgs_params.mem_size = 256;
            lbfgs_params.past = 3;
            lbfgs_params.min_step = 1.0e-32;
            lbfgs_params.max_step = 100;
            lbfgs_params.g_epsilon = 1e-6;
            lbfgs_params.delta = relCostTol;
            lbfgs_params.max_iterations = 100000;

            lbfgs_params.use_alm = true;
            lbfgs_params.alm_delta = 1.;
            lbfgs_params.mu_rho = &mu_rho;
            lbfgs_params.projection_on_symm_cone = &projection_on_symm_cone;
            lbfgs_params.equal_criterion = &equal_criterion;
            lbfgs_params.inequal_criterion = &inequal_criterion;
            lbfgs_params.conic_crit = DBL_MAX;
            lbfgs_params.equal_crit = DBL_MAX;
            lbfgs_params.inequal_crit = DBL_MAX;

            rho = 1;
            mu_rho = mu / rho;
            lambda_rho = lambda / rho;
            sigma_rho = sigma / rho;

            double minCost;
            const int MAX_COUNT = 100000;
            int count = 0;
            while (++count < MAX_COUNT)
            {
                int ret = lbfgs::lbfgs_optimize(x,
                                                minCost,
                                                &Topp::costFunction,
                                                nullptr,
                                                nullptr,
                                                this,
                                                lbfgs_params);
                if (ret < 0)
                {
                    std::cout << "\033[32m" << lbfgs::lbfgs_strerror(ret) << "\033[0m" << std::endl;
                    return ret;
                }
                // std::cout << "\033[32m"
                //           << count << " alm iteration, res: " << ret << "\033[0m" << std::endl;
                if (ret >= 0 && lbfgs_params.conic_crit < relCostTol && lbfgs_params.equal_crit < relCostTol && lbfgs_params.inequal_crit < relCostTol)
                {
                    bs = x.segment(0, K + 1);
                    as = x.segment(K + 1, K);
                    return ret;
                }
                mu = projection_on_symm_cone * rho;
                lambda += equal_criterion * rho;
                rho = std::min((1. + gamma_) * rho, beta_);
                mu_rho = mu / rho;
                lambda_rho = lambda / rho;
                sigma_rho = sigma / rho;
                lbfgs_params.alm_delta = lbfgs_params.alm_delta / 2;
            }

            if (count == MAX_COUNT)
            {
                std::cout << "Failed: topp ALM count max" << std::endl;
                return -1;
            }
            return minCost;
        }
    };
}

#endif
