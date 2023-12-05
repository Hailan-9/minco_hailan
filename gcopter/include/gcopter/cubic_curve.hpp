#ifndef CUBIC_CURVE_HPP
#define CUBIC_CURVE_HPP

#include <Eigen/Eigen>

#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
// 三次多项式
// 三次多项式
// 三次多项式

// 多项式的时间使用相对时间，每段的时间都是从0开始
class CubicPolynomial
{
private:
    double duration;
    // 两个维度x、y
    Eigen::Matrix<double, 2, 4> coeffMat;

public:
    CubicPolynomial() = default;

    CubicPolynomial(double dur, const Eigen::Matrix<double, 2, 4> &cMat)
        : duration(dur), coeffMat(cMat) {}

    inline int getDim() const
    {
        return 2;
    }

    inline int getDegree() const
    {
        return 3;
    }

    inline double getDuration() const
    {
        return duration;
    }

    inline const Eigen::Matrix<double, 2, 4> &getCoeffMat() const
    {
        return coeffMat;
    }
    // x、y两个维度的位置 系数对应的阶次从高到低
    inline Eigen::Vector2d getPos(const double &t) const
    {
        return coeffMat.col(3) + t * (coeffMat.col(2) + t * (coeffMat.col(1) + t * coeffMat.col(0)));
    }
    inline Eigen::Vector2d getVel(const double &t) const
    {
        return coeffMat.col(2) + t * 2 * coeffMat.col(1) + t * t * 3 * coeffMat.col(0);
    }
    inline Eigen::Vector2d getAcc(const double &t) const
    {
        return 2 * coeffMat.col(1) + t * 6 * coeffMat.col(0);
    }
};
// 三次曲线
class CubicCurve
{
private:
    // 多段三次多项式曲线
    typedef std::vector<CubicPolynomial> Pieces;
    Pieces pieces;

public:
    CubicCurve() = default;

    CubicCurve(const std::vector<double> &durs,
               const std::vector<Eigen::Matrix<double, 2, 4>> &cMats)
    {
        const int N = std::min(durs.size(), cMats.size());
        pieces.reserve(N);
        for (int i = 0; i < N; ++i)
        {
            pieces.emplace_back(durs[i], cMats[i]);
        }
    }
    // 返回三次多项式的段数，即轨迹段数
    inline int getPieceNum() const
    {
        return pieces.size();
    }

    inline Eigen::VectorXd getDurations() const
    {
        const int N = getPieceNum();
        Eigen::VectorXd durations(N);
        for (int i = 0; i < N; ++i)
        {
            durations(i) = pieces[i].getDuration();
        }
        return durations;
    }

    inline double getTotalDuration() const
    {
        const int N = getPieceNum();
        double totalDuration = 0.0;
        for (int i = 0; i < N; ++i)
        {
            totalDuration += pieces[i].getDuration();
        }
        return totalDuration;
    }
    // 得到航迹中间点，即每段轨迹的对应的两个点（2d）
    inline Eigen::Matrix2Xd getPositions() const
    {
        const int N = getPieceNum();
        Eigen::Matrix2Xd positions(2, N + 1);
        for (int i = 0; i < N; ++i)
        {
            positions.col(i) = pieces[i].getCoeffMat().col(3);
        }
        positions.col(N) = pieces[N - 1].getPos(pieces[N - 1].getDuration());
        return positions;
    }
    // 成员函数实现运算符重载
    inline const CubicPolynomial &operator[](int i) const
    {
        return pieces[i];
    }

    inline CubicPolynomial &operator[](int i)
    {
        return pieces[i];
    }

    inline void clear(void)
    {
        pieces.clear();
        return;
    }

    inline Pieces::const_iterator begin() const
    {
        return pieces.begin();
    }

    inline Pieces::const_iterator end() const
    {
        return pieces.end();
    }

    inline Pieces::iterator begin()
    {
        return pieces.begin();
    }

    inline Pieces::iterator end()
    {
        return pieces.end();
    }

    inline void reserve(const int &n)
    {
        pieces.reserve(n);
        return;
    }

    inline void emplace_back(const CubicPolynomial &piece)
    {
        pieces.emplace_back(piece);
        return;
    }

    inline void emplace_back(const double &dur,
                             const Eigen::Matrix<double, 2, 4> &cMat)
    {
        pieces.emplace_back(dur, cMat);
        return;
    }

    inline void append(const CubicCurve &traj)
    {
        pieces.insert(pieces.end(), traj.begin(), traj.end());
        return;
    }
    // 定位一下传入的时间参数在哪段三次多项式上面
    inline int locatePieceIdx(double &t) const
    {
        const int N = getPieceNum();
        int idx;
        double dur;
        // 先赋值后比较
        for (idx = 0;
             idx < N &&
             t > (dur = pieces[idx].getDuration());
             idx++)
        {
            t -= dur;
        }
        // 会进入这一项吗？？？
        if (idx == N)
        {
            idx--;
            t += pieces[idx].getDuration();
        }
        return idx;
    }

    inline Eigen::Vector2d getPos(double t) const
    {
        const int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getPos(t);
    }

    inline void getPosVelAcc(double t, Eigen::Vector2d &pos, Eigen::Vector2d &vel, Eigen::Vector2d &acc) const
    {
        const int pieceIdx = locatePieceIdx(t);
        pos = pieces[pieceIdx].getPos(t);
        vel = pieces[pieceIdx].getVel(t);
        acc = pieces[pieceIdx].getAcc(t);
    }

    inline Eigen::Vector2d getJuncPos(const int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getCoeffMat().col(3);
        }
        else
        {
            return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
        }
    }
};

#endif
