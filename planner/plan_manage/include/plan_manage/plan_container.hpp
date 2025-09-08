#ifndef _PLAN_CONTAINER_H_
#define _PLAN_CONTAINER_H_

#include <Eigen/Eigen>
#include <vector>
#include <bspline_opt/uniform_bspline.h>

using std::vector;

namespace ego_planner
{
  struct PlanParameters
  {
	//规划过程相关参数
    double max_vel_, max_acc_, max_jerk_; // 物理约束、最大速度、最大加速度、最大加加速度（Jerk_)
    double ctrl_pt_dist;                  // B样条曲线中的控制点间距
    double feasibility_tolerance_;        // 允许超过物理约束的比例
    double planning_horizen_; 			// 规划视野范围

    /* 各个阶段的耗时 */
    double time_search_ = 0.0;
    double time_optimize_ = 0.0;
    double time_adjust_ = 0.0;
  };

  struct LocalTrajData
  {
    /* 局部轨迹数据 */

    int traj_id_; // 轨迹ID
    double duration_; // 轨迹持续时间
    double global_time_offset; // 局部轨迹结束时，时间无法匹配全局轨迹，需要偏移
    Eigen::Vector3d start_pos_; // 开始位置
    UniformBspline position_traj_, velocity_traj_, acceleration_traj_;
  };

} // namespace ego_planner

#endif
