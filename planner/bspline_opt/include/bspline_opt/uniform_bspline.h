#ifndef _UNIFORM_BSPLINE_H_
#define _UNIFORM_BSPLINE_H_

#include <Eigen/Eigen> // 包含Eigen库，用于矩阵和向量操作
#include <algorithm> // 包含C++标准库中的算法
#include <iostream> // 包含C++标准库中的输入输出流

using namespace std; // 使用标准命名空间

namespace ego_planner
{
  // 非均匀B样条曲线的实现，支持不同维度
  // 它也表示均匀B样条曲线，这是一种特殊的非均匀B样条曲线
  class UniformBspline
  {
  private:
    // 不同维度的B样条曲线的控制点矩阵
    // 每一行代表一个单独的控制点
    // 维度由列数决定
    // 例如：3D空间中的N个点 -> Nx3矩阵
    Eigen::MatrixXd control_points_;

    int p_; // B样条曲线的阶数（degree）
    int n_; // 控制点的数量（n+1）
    int m_;  
    Eigen::VectorXd u_; // 节点向量（knots vector）
    double interval_;   // 节点间隔（knot span \delta t）

    // 计算B样条曲线的导数控制点
    Eigen::MatrixXd getDerivativeControlPoints();

    double limit_vel_; // 速度限制
    double limit_acc_; // 加速度限制
    double limit_ratio_; // 时间调整比例
    double feasibility_tolerance_; // 可行性容差

  public:
    // 默认构造函数
    UniformBspline() {}

    // 带参数的构造函数，初始化B样条曲线
    UniformBspline(const Eigen::MatrixXd &points, const int &order, const double &interval);

    // 析构函数
    ~UniformBspline();

    // 获取控制点矩阵
    Eigen::MatrixXd get_control_points(void) { return control_points_; }

    // 初始化为均匀B样条曲线
    void setUniformBspline(const Eigen::MatrixXd &points, const int &order, const double &interval);

    // 获取或设置基本的B样条信息

    // 设置结点向量
    void setKnot(const Eigen::VectorXd &knot);

    // 获取结点向量
    Eigen::VectorXd getKnot();

    // 获取控制点矩阵
    Eigen::MatrixXd getControlPoint();

    // 获取结点间隔
    double getInterval();

    // 获取时间范围
    bool getTimeSpan(double &um, double &um_p);

    // 使用De Boor算法计算位置 / 导数

    // 使用结点参数u计算B样条曲线的位置，u \in [up, u_mp]
    Eigen::VectorXd evaluateDeBoor(const double &u);

    // 使用时间参数t计算B样条曲线的位置，t \in [0, duration] 计算开始后t秒的点
    inline Eigen::VectorXd evaluateDeBoorT(const double &t) { return evaluateDeBoor(t + u_(p_)); }

    // 获取当前B样条曲线的导数
    UniformBspline getDerivative();

    // 3D B样条曲线插值，对point_set中的点进行插值，并考虑边界速度和加速度约束
    // 输入：(K+2)个点以及边界速度/加速度；时间间隔ts
    // 输出：(K+6)个控制点ctrl_pts
    static void parameterizeToBspline(const double &ts, const vector<Eigen::Vector3d> &point_set,
                                      const vector<Eigen::Vector3d> &start_end_derivative,
                                      Eigen::MatrixXd &ctrl_pts);

    /* 检查可行性并调整时间 */

    // 设置物理限制，包括速度、加速度和容差
    void setPhysicalLimits(const double &vel, const double &acc, const double &tolerance);

    // 检查B样条曲线的可行性，并返回调整比例ratio
    bool checkFeasibility(double &ratio, bool show = false);

    // 根据给定的比例ratio延长B样条曲线的时间
    void lengthenTime(const double &ratio);

    /* 用于性能评估 */

    // 获取B样条曲线的总时间
    double getTimeSum();

    // 获取B样条曲线的长度
    double getLength(const double &res = 0.01);

    // 获取B样条曲线的 jerk（加加速度）
    double getJerk();

    // 获取B样条曲线的平均速度和最大速度
    void getMeanAndMaxVel(double &mean_v, double &max_v);

    // 获取B样条曲线的平均加速度和最大加速度
    void getMeanAndMaxAcc(double &mean_a, double &max_a);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // 确保Eigen对象正确对齐
  };
} // namespace ego_planner
#endif
