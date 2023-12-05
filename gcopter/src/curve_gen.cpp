#include "misc/visualizer.hpp"
#include "gcopter/cubic_curve.hpp"
#include "gcopter/cubic_spline.hpp"
#include "gcopter/path_smoother.hpp"
#include "gcopter/topp.hpp"

#include <ros/ros.h>
#include <ros/console.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <random>
#define NUMMOD 10000
#define NUMDEV 10000.0
struct Config
{
    std::string targetTopic;
    double penaltyWeight;
    Eigen::Matrix3Xd circleObs;
    double pieceLength;
    double relCostTol;
    std::vector<Eigen::MatrixX4d> hPolys;
    // 障碍物
    std::vector<Eigen::MatrixX3d> hPolys2d;
    double v_max, a_max;
    double v_start, v_end;
    int sub_pieces_num;

    Config(const ros::NodeHandle &nh_priv)
    {
        std::vector<double> circleObsVec;
        int xy_bound = 25, uv_bound = 5, planes_num = 10, polys_num = 20;
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("PenaltyWeight", penaltyWeight);
        nh_priv.getParam("CircleObs", circleObsVec);
        nh_priv.getParam("PieceLength", pieceLength);
        nh_priv.getParam("RelCostTol", relCostTol);
        nh_priv.getParam("PolysNum", polys_num);
        nh_priv.getParam("PolyPlaneNum", planes_num);
        nh_priv.getParam("PolyXYBound", xy_bound);
        nh_priv.getParam("PolyUVBound", uv_bound);
        nh_priv.getParam("v_max", v_max);
        nh_priv.getParam("a_max", a_max);
        nh_priv.getParam("v_start", v_start);
        nh_priv.getParam("v_end", v_end);
        nh_priv.getParam("sub_pieces_num", sub_pieces_num);
        // circleObs = Eigen::Map<const Eigen::Matrix<double, 3, -1, Eigen::ColMajor>>(
        //     circleObsVec.data(), 3, circleObsVec.size() / 3);

        /**generate hPolys obs
         * center (x0,y0), a plane vector from center is  (u,v), which are random
         * plane: [u, v, -(u^2 + v^2 + ux0 + uy0)] * [x, y, 1]^T <= 0
         **/
        double x0, y0, u, v;
        Eigen::Matrix<double, 2, 4> top_below_plane; // for 3d display
        top_below_plane.setZero();
        top_below_plane(0, 2) = 1;
        top_below_plane(0, 3) = -1;
        top_below_plane(1, 2) = -1;
        top_below_plane(1, 3) = -1;

        srand((unsigned)std::time(NULL));
        circleObs.resize(3, polys_num);
        for (int j = 0; j < polys_num; ++j)
        {
            Eigen::MatrixX4d hPoly3d;
            Eigen::MatrixX3d hPoly2d;
            hPoly3d.resize(planes_num + 2, 4);
            hPoly2d.resize(planes_num, 3);
            x0 = rand() % (2 * xy_bound * NUMMOD) / NUMDEV - xy_bound;
            y0 = rand() % (2 * xy_bound * NUMMOD) / NUMDEV - xy_bound;
            circleObs(0, j) = x0;
            circleObs(1, j) = y0;
            circleObs(2, j) = 0.5;
            for (int i = 0; i < planes_num; i++)
            {
                u = rand() % (2 * uv_bound * NUMMOD) / NUMDEV - uv_bound;
                v = rand() % (2 * uv_bound * NUMMOD) / NUMDEV - uv_bound;
                hPoly2d(i, 0) = u;
                hPoly2d(i, 1) = v;
                hPoly2d(i, 2) = -(u * u + v * v + u * x0 + v * y0);
                // std::cout << "x0 = " << x0 << ", y0 = " << y0 << std::endl;
                // std::cout << "u = " << u << ", v = " << v << std::endl;
            }
            hPoly3d.setZero();
            hPoly3d.block(0, 0, planes_num, 2) = hPoly2d.block(0, 0, planes_num, 2);
            hPoly3d.block(0, 3, planes_num, 1) = hPoly2d.block(0, 2, planes_num, 1);
            hPoly3d.block(planes_num, 0, 2, 4) = top_below_plane;
            hPolys.push_back(hPoly3d);
            hPolys2d.push_back(hPoly2d);
        }
        // Eigen::Matrix<double, -1, 2> A = hPolys2d.at(0).leftCols(2);
        // Eigen::Matrix<double, -1, 1> b = hPolys2d.at(0).rightCols(1);
        // Eigen::Matrix<double, 2, 1> x(x0, y0);
        // std::cout << "cons precision: " << (A * x - b).maxCoeff() << std::endl;
    }
};

class CurveGen
{
private:
    Config config;
    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;
    std::vector<Eigen::Vector2d> startGoal;

    CubicCurve curve;

public:
    CurveGen(ros::NodeHandle &nh_)
        : config(ros::NodeHandle("~")),
          nh(nh_),
          visualizer(nh)
    {
        targetSub = nh.subscribe(config.targetTopic, 1, &CurveGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    inline void vizObs()
    {
        // visualizer.visualizeDisks(config.circleObs);
        visualizer.visualizePolytope(config.hPolys);
    }

    // 主要看这个函数，奥利给！！！
    inline void plan()
    {
        if (startGoal.size() == 2)
        {
            // 计算轨迹段数
            const int N = (startGoal.back() - startGoal.front()).norm() / config.pieceLength;
            Eigen::Matrix2Xd innerPoints(2, N - 1); // 轨迹内部点，即去除起点终点俩端点后，剩余的航迹中间点
            for (int i = 0; i < N - 1; ++i) // 中间内部点，均匀分布
            {
                innerPoints.col(i) = (startGoal.back() - startGoal.front()) * (i + 1.0) / N + startGoal.front();
            }

            // 1.路径平滑
            path_smoother::PathSmoother pathSmoother;
            // pathSmoother.setup(startGoal.front(), startGoal.back(), N, config.circleObs, config.penaltyWeight);
            pathSmoother.setup(startGoal.front(), startGoal.back(), N, config.hPolys2d, config.penaltyWeight);

            CubicCurve curve;
            // 代价容许 tolerance
            if (std::isinf(pathSmoother.optimize(curve, innerPoints, config.relCostTol)))
            {
                return;
            }
            if (curve.getPieceNum() > 0)
            {
                visualizer.visualize(curve);
            }

            /***  topp  ***/
            const int RESO = config.sub_pieces_num;
            topp::Topp topp_solver(curve, RESO);
            Eigen::VectorXd bs(RESO * N + 1);
            Eigen::VectorXd as(RESO * N);

            int ret = topp_solver.optimize(bs, as, config.v_max, config.a_max, config.v_start, config.v_end, config.relCostTol);
            // topp_solver.printInfo();

            if (ret >= 0)
            {
                std::cout << "[curv_gen] TOPP success" << std::endl;
                std::cout << "bs:\n"
                          << bs;
                std::cout << "\n as: \n"
                          << as << std::endl;

                // pub result
                visualizer.pubVectorXd(bs, "bs");
                visualizer.pubVectorXd(as, "as");
                visualizer.pubVectorXd(topp_solver.s, "s");
                visualizer.pubVectorXd(topp_solver.ax, "accx");
                visualizer.pubVectorXd(topp_solver.ay, "accy");
            }
        }
    }

    inline void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {

        if (startGoal.size() >= 2)
        {
            startGoal.clear();
        }

        startGoal.emplace_back(msg->pose.position.x, msg->pose.position.y);

        plan();

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "curve_gen_node");
    ros::NodeHandle nh_;

    CurveGen curveGen(nh_);

    ros::Duration(2.0).sleep();

    curveGen.vizObs();

    ros::spin();

    return 0;
}
