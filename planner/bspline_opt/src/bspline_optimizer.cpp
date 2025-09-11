#include "bspline_opt/bspline_optimizer.h"
#include "bspline_opt/gradient_descent_optimizer.h"
// using namespace std;

namespace ego_planner
{

  void BsplineOptimizer::setParam()
  {
    lambda1_ = 10.0;
    lambda2_ = 0.5;
    lambda3_ = 0.1;
    lambda4_ = 1.0;
    dist0_   = 0.5;   // 安全距离
    max_vel_  = 1.0;  // 最大速度
    max_acc_  = 1.0;  // 最大加速度
    order_    = 3;    // b样条阶数

  }

  void BsplineOptimizer::setEnvironment(const shared_ptr<SDFMap> &env)
  {
    this->grid_map_ = env;
  }

  void BsplineOptimizer::setControlPoints(const Eigen::MatrixXd &points)
  {
    cps_.points = points;
  }

  void BsplineOptimizer::setBsplineInterval(const double &ts) { bspline_interval_ = ts; }

  /* This function is very similar to check_collision_and_rebound(). 
   * It was written separately, just because I did it once and it has been running stably since March 2020.
   * But I will merge then someday.*/
   
  //初始化控制点(传入B样条的控制点）
  std::vector<std::vector<Eigen::Vector3d>> BsplineOptimizer::initControlPoints(Eigen::MatrixXd &init_points, bool flag_first_init /*= true*/)
  {

    if (flag_first_init)
    {
      cps_.clearance = dist0_;
      cps_.resize(init_points.cols());
      cps_.points = init_points; //传入控制点
    }

    /*** 根据障碍物给初始轨迹分段 ***/
    constexpr int ENOUGH_INTERVAL = 2;
    // 计算步长 分辨率 / (首尾路径点距离/路径点数量) / 2
    double step_size = grid_map_->getResolution() / ((init_points.col(0) - init_points.rightCols(1)).norm() / (init_points.cols() - 1)) / 2;
    int in_id, out_id;
    vector<std::pair<int, int>> segment_ids; //存储分段索引
    int same_occ_state_times = ENOUGH_INTERVAL + 1; // 检测是不是同一个障碍物的状态计数器
    bool occ, last_occ = false; //是否是障碍物点，上次是否是障碍物点
    bool flag_got_start = false, flag_got_end = false, flag_got_end_maybe = false; //标志位，是否找到起点，终点，终点可能
    int i_end = (int)init_points.cols() - order_ - ((int)init_points.cols() - 2 * order_) / 3; // 只检查前面约 2/3 的点 
    for (int i = order_; i <= i_end; ++i) //遍历
    {
      for (double a = 1.0; a >= 0.0; a -= step_size) //
      {
        // 插值迭代控制点连线
        occ = grid_map_->getInflateOccupancy(a * init_points.col(i - 1) + (1 - a) * init_points.col(i));
        // cout << setprecision(5);
        // cout << (a * init_points.col(i-1) + (1-a) * init_points.col(i)).transpose() << " occ1=" << occ << endl;

        //查找分段的头部和尾部
        if (occ && !last_occ) // 进入障碍物 
        {
          if (same_occ_state_times > ENOUGH_INTERVAL || i == order_)//开始位置或者具备足够的间隔
          {
            in_id = i - 1;//标记分段开始索引
            flag_got_start = true;
          }
          same_occ_state_times = 0; //重制计数器
          flag_got_end_maybe = false; // terminate in advance
        }
        else if (!occ && last_occ) // 离开障碍物
        {
          out_id = i; //标记分段结束索引
          flag_got_end_maybe = true;
          same_occ_state_times = 0;
        }
        else //连续同障碍物索引状态 （连续为空闲 或连续为占用）
        {
          ++same_occ_state_times;
        }

        // 可能结尾状态时检查到了最后的控制点
        if (flag_got_end_maybe && (same_occ_state_times > ENOUGH_INTERVAL || (i == (int)init_points.cols() - order_))) // 这部分不会触发，因为i_end 已经限制了检查范围
        {
          flag_got_end_maybe = false;
          flag_got_end = true;
        }

        last_occ = occ; //更新上次状态
        //找到段尾和段头，记录当前分段
        if (flag_got_start && flag_got_end)
        {
          flag_got_start = false;
          flag_got_end = false;
          segment_ids.push_back(std::pair<int, int>(in_id, out_id));
        }
      }
    }

    /*** a star search ***/
    // 对找到障碍物的分段执行A*搜索获得局部路径
    vector<vector<Eigen::Vector3d>> a_star_pathes;
    for (size_t i = 0; i < segment_ids.size(); ++i)
    {
      //cout << "in=" << in.transpose() << " out=" << out.transpose() << endl;
      Eigen::Vector3d in(init_points.col(segment_ids[i].first)), out(init_points.col(segment_ids[i].second));
      if (a_star_->AstarSearch(/*(in-out).norm()/10+0.05*/ 0.1, in, out))
      {
        a_star_pathes.push_back(a_star_->getPath());
      }
      else
      {
        printf("a star error, force return!\n");
        return a_star_pathes;
      }
    }

    /*** 控制点边界计算 ***/
    int id_low_bound, id_up_bound;
    vector<std::pair<int, int>> bounds(segment_ids.size());
    for (size_t i = 0; i < segment_ids.size(); i++)
    {

      if (i == 0) // 第一个分段
      {
        id_low_bound = order_; // 定为首个控制点
        if (segment_ids.size() > 1) // 不止一段
        {
          // 将边界设置在当前分段结尾和下一分段开始之间间隔的中间控制点
          // +-1.0f 用于修正部分计算问题？
          id_up_bound = (int)(((segment_ids[0].second + segment_ids[1].first) - 1.0f) / 2); // id_up_bound : -1.0f fix()
        }
        else
        {
          // 只有一段 设置为最后一个控制点
          id_up_bound = init_points.cols() - order_ - 1;
        }
      }
      else if (i == segment_ids.size() - 1) // 最后一个分段
      {
        id_low_bound = (int)(((segment_ids[i].first + segment_ids[i - 1].second) + 1.0f) / 2); // 与上一段结尾点间隔的中间点
        id_up_bound = init_points.cols() - order_ - 1; //设为最后一个控制点
      }
      else
      {
        id_low_bound = (int)(((segment_ids[i].first + segment_ids[i - 1].second) + 1.0f) / 2); // id_low_bound : +1.0f ceil()
        id_up_bound = (int)(((segment_ids[i].second + segment_ids[i + 1].first) - 1.0f) / 2);  // id_up_bound : -1.0f fix()
      }

      bounds[i] = std::pair<int, int>(id_low_bound, id_up_bound); //存入bounds
    }

    // cout << "+++++++++" << endl;
    // for ( int j=0; j<bounds.size(); ++j )
    // {
    //   cout << bounds[j].first << "  " << bounds[j].second << endl;
    // }

    /*** Adjust segment length ***/
    //按照点的数量调整分段长度 由于 MINIMUM_PERCENT = 0.0;这部分没有效果
    vector<std::pair<int, int>> final_segment_ids(segment_ids.size());
    constexpr double MINIMUM_PERCENT = 0.0; // Each segment is guaranteed to have sufficient points to generate sufficient thrust 保证每个部分都有足够的点来产生足够的推力 最小点数百分比 此处没有用
    int minimum_points = round(init_points.cols() * MINIMUM_PERCENT), num_points; //最小点数量
    for (size_t i = 0; i < segment_ids.size(); i++)
    {
      /*** Adjust segment length ***/
      num_points = segment_ids[i].second - segment_ids[i].first + 1; //控制点数量
      //cout << "i = " << i << " first = " << segment_ids[i].first << " second = " << segment_ids[i].second << endl;
      if (num_points < minimum_points)
      {
        double add_points_each_side = (int)(((minimum_points - num_points) + 1.0f) / 2);

        final_segment_ids[i].first = segment_ids[i].first - add_points_each_side >= bounds[i].first ? segment_ids[i].first - add_points_each_side : bounds[i].first;

        final_segment_ids[i].second = segment_ids[i].second + add_points_each_side <= bounds[i].second ? segment_ids[i].second + add_points_each_side : bounds[i].second;
      }
      else //复制到final_segment_ids
      {
        final_segment_ids[i].first = segment_ids[i].first;
        final_segment_ids[i].second = segment_ids[i].second;
      }

      //cout << "final:" << "i = " << i << " first = " << final_segment_ids[i].first << " second = " << final_segment_ids[i].second << endl;
    }

    /*** Assign data to each segment ***/
    // 给每个分段分配数据
    for (size_t i = 0; i < segment_ids.size(); i++)//遍历每个分段
    {
      // 第一步 将分段内控制点初始化为未标记
      for (int j = final_segment_ids[i].first; j <= final_segment_ids[i].second; ++j)
        cps_.flag_temp[j] = false; // 标记为未使用

      int got_intersection_id = -1; //j 对应的插值点索引
      for (int j = segment_ids[i].first + 1; j < segment_ids[i].second; ++j) //遍历除了端点以外的控制点 （这里没有处理没有中间控制点的情况，在后续处理）
      {
        // 控制点运动方向
        Eigen::Vector3d ctrl_pts_law(cps_.points.col(j + 1) - cps_.points.col(j - 1)), intersection_point;
        //Astart_id 标记为中间点的索引
        //注释为选择最远的A*点会更好，但是需要更多的计算量
        int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id; // Let "Astar_id = id_of_the_most_far_away_Astar_point" will be better, but it needs more computation
        //叉积 用于在符号变化时检测交点
        double val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot(ctrl_pts_law), last_val = val; //A*点相对于ctrl_pts_law方向的投影值
        while (Astar_id >= 0 && Astar_id < (int)a_star_pathes[i].size()) //迭代A*路径
        {
          last_Astar_id = Astar_id;

          //沿着
          if (val >= 0)
            --Astar_id;
          else
            ++Astar_id;

          val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot(ctrl_pts_law);

          // 当投影方向发生变化时，即默认此时A*上点Astar_id点 从进程上对应着原控制点
          if (val * last_val <= 0 && (abs(val) > 0 || abs(last_val) > 0)) // val = last_val = 0.0 is not allowed 排除垂直夹角的情况
          {
            //使用线形插值计算插值点
            intersection_point =
                a_star_pathes[i][Astar_id] +
                ((a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id]) * //A*段方向
                  (
                    ctrl_pts_law.dot(cps_.points.col(j) - a_star_pathes[i][Astar_id]) //相邻控制点的运动向量在(接近垂直方向上的)投影长度
                    / 
                    ctrl_pts_law.dot(a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id]) // A*段方向的投影长度
                  ) // = t
                );

            //cout << "i=" << i << " j=" << j << " Astar_id=" << Astar_id << " last_Astar_id=" << last_Astar_id << " intersection_point = " << intersection_point.transpose() << endl;

            got_intersection_id = j; //记录插值点对应的控制点索引
            break;
          }
        }

        //如果存在插值点
        if (got_intersection_id >= 0)
        {
          cps_.flag_temp[j] = true; //对应的控制点标记为已修改
          double length = (intersection_point - cps_.points.col(j)).norm(); //插值点到原控制点的距离
          if (length > 1e-5) //排除过近的情况
          {
            for (double a = length; a >= 0.0; a -= grid_map_->getResolution())
            {
              occ = grid_map_->getInflateOccupancy((a / length) * intersection_point + (1 - a / length) * cps_.points.col(j)); //在原点和插值点连线上线形插值并判断障碍物情况 //由插值点向原点方向插值

              if (occ || a < grid_map_->getResolution()) // 障碍物或过近
              {
                if (occ)
                  a += grid_map_->getResolution(); //取上一次的位置
                //将确定的点保存
                cps_.base_point[j].push_back((a / length) * intersection_point + (1 - a / length) * cps_.points.col(j));
                cps_.direction[j].push_back((intersection_point - cps_.points.col(j)).normalized());
                break;
              }
            }
          }
        }
      }

      /* Corner case: the segment length is too short. Here the control points may outside the A* path, leading to opposite gradient direction. So I have to take special care of it */
      //如果该分段没有中间控制点
      //取一个中点middle_point进行构造和插值（计算大致与上面相同）
      if (segment_ids[i].second - segment_ids[i].first == 1)
      {
        Eigen::Vector3d ctrl_pts_law(cps_.points.col(segment_ids[i].second) - cps_.points.col(segment_ids[i].first)), intersection_point;//计算运动方向
        Eigen::Vector3d middle_point = (cps_.points.col(segment_ids[i].second) + cps_.points.col(segment_ids[i].first)) / 2; //取中间点
        int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id; // Let "Astar_id = id_of_the_most_far_away_Astar_point" will be better, but it needs more computation
        double val = (a_star_pathes[i][Astar_id] - middle_point).dot(ctrl_pts_law), last_val = val;
        while (Astar_id >= 0 && Astar_id < (int)a_star_pathes[i].size())
        {
          last_Astar_id = Astar_id;

          if (val >= 0)
            --Astar_id;
          else
            ++Astar_id;

          val = (a_star_pathes[i][Astar_id] - middle_point).dot(ctrl_pts_law);

          if (val * last_val <= 0 && (abs(val) > 0 || abs(last_val) > 0)) // val = last_val = 0.0 is not allowed
          {
            intersection_point =
                a_star_pathes[i][Astar_id] +
                ((a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id]) *
                 (ctrl_pts_law.dot(middle_point - a_star_pathes[i][Astar_id]) / ctrl_pts_law.dot(a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id])) // = t
                );

            if ((intersection_point - middle_point).norm() > 0.01) // 1cm.
            {
              cps_.flag_temp[segment_ids[i].first] = true;
              cps_.base_point[segment_ids[i].first].push_back(cps_.points.col(segment_ids[i].first));
              cps_.direction[segment_ids[i].first].push_back((intersection_point - middle_point).normalized());

              got_intersection_id = segment_ids[i].first; //插值点id标记为分段的第一个点
            }
            break;
          }
        }
      }

      //step 3
      if (got_intersection_id >= 0)
      {
        for (int j = got_intersection_id + 1; j <= final_segment_ids[i].second; ++j)
          if (!cps_.flag_temp[j])
          {
            cps_.base_point[j].push_back(cps_.base_point[j - 1].back());
            cps_.direction[j].push_back(cps_.direction[j - 1].back());
          }

        for (int j = got_intersection_id - 1; j >= final_segment_ids[i].first; --j)
          if (!cps_.flag_temp[j])
          {
            cps_.base_point[j].push_back(cps_.base_point[j + 1].back());
            cps_.direction[j].push_back(cps_.direction[j + 1].back());
          }
      }
      else
      {
        // Just ignore, it does not matter ^_^.
        // ROS_ERROR("Failed to generate direction! segment_id=%d", i);
      }
    }

    return a_star_pathes;
  }

  int BsplineOptimizer::earlyExit(void *func_data, const double *x, const double *g, const double fx, const double xnorm, const double gnorm, const double step, int n, int k, int ls)
  {
    BsplineOptimizer *opt = reinterpret_cast<BsplineOptimizer *>(func_data);
    // cout << "k=" << k << endl;
    // cout << "opt->flag_continue_to_optimize_=" << opt->flag_continue_to_optimize_ << endl;
    return (opt->force_stop_type_ == STOP_FOR_ERROR || opt->force_stop_type_ == STOP_FOR_REBOUND);
  }

  // 包装函数，用于传入lbfgs
  double BsplineOptimizer::costFunctionRebound(void *func_data, const double *x, double *grad, const int n)
  {
    BsplineOptimizer *opt = reinterpret_cast<BsplineOptimizer *>(func_data);

    double cost;
    opt->combineCostRebound(x, grad, cost, n);

    opt->iter_num_ += 1;
    return cost;
  }

  // 包装函数，用于传入lbfgs
  double BsplineOptimizer::costFunctionRefine(void *func_data, const double *x, double *grad, const int n)
  {
    BsplineOptimizer *opt = reinterpret_cast<BsplineOptimizer *>(func_data);

    double cost;
    opt->combineCostRefine(x, grad, cost, n);

    opt->iter_num_ += 1;
    return cost;
  }

  //计算 障碍物距离代价和梯度
  void BsplineOptimizer::calcDistanceCostRebound(const Eigen::MatrixXd &q, double &cost,
                                                 Eigen::MatrixXd &gradient, int iter_num, double smoothness_cost)
  {
    cost = 0.0;
    int end_idx = q.cols() - order_;
    double demarcation = cps_.clearance; //安全距离
    double  a = 3 * demarcation, 
            b = -3 * pow(demarcation, 2), 
            c = pow(demarcation, 3);

    force_stop_type_ = DONT_STOP;
    //确保足够平滑？
    if (iter_num > 3 && smoothness_cost / (cps_.size - 2 * order_) < 0.1) // 0.1 is an experimental value that indicates the trajectory is smooth enough.
    {
      check_collision_and_rebound();
    }

    /*** calculate distance cost and gradient ***/
    for (auto i = order_; i < end_idx; ++i)
    {
      for (size_t j = 0; j < cps_.direction[i].size(); ++j)
      {
        double dist = (cps_.points.col(i) - cps_.base_point[i][j]).dot(cps_.direction[i][j]);
        double dist_err = cps_.clearance - dist;
        Eigen::Vector3d dist_grad = cps_.direction[i][j];

        if (dist_err < 0)
        {
          /* do nothing */
        }
        else if (dist_err < demarcation)
        {
          cost += pow(dist_err, 3);
          gradient.col(i) += -3.0 * dist_err * dist_err * dist_grad;
        }
        else
        {
          cost += a * dist_err * dist_err + b * dist_err + c;
          gradient.col(i) += -(2.0 * a * dist_err + b) * dist_grad;
        }
      }
    }
  }

  void BsplineOptimizer::calcFitnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)
  {

    cost = 0.0;

    int end_idx = q.cols() - order_;

    // def: f = |x*v|^2/a^2 + |x×v|^2/b^2
    double a2 = 25, b2 = 1;
    for (auto i = order_ - 1; i < end_idx + 1; ++i)
    {
      Eigen::Vector3d x = (q.col(i - 1) + 4 * q.col(i) + q.col(i + 1)) / 6.0 - ref_pts_[i - 1];
      Eigen::Vector3d v = (ref_pts_[i] - ref_pts_[i - 2]).normalized();

      double xdotv = x.dot(v);
      Eigen::Vector3d xcrossv = x.cross(v);

      double f = pow((xdotv), 2) / a2 + pow(xcrossv.norm(), 2) / b2;
      cost += f;

      Eigen::Matrix3d m;
      m << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
      Eigen::Vector3d df_dx = 2 * xdotv / a2 * v + 2 / b2 * m * xcrossv;

      gradient.col(i - 1) += df_dx / 6;
      gradient.col(i) += 4 * df_dx / 6;
      gradient.col(i + 1) += df_dx / 6;
    }
  }

  //计算平滑代价和梯度
  void BsplineOptimizer::calcSmoothnessCost(const Eigen::MatrixXd &q, double &cost,
                                            Eigen::MatrixXd &gradient, bool falg_use_jerk /* = true*/)
  {

    cost = 0.0;

    if (falg_use_jerk) //使用加加速度
    {
      Eigen::Vector3d jerk, temp_j;

      for (int i = 0; i < q.cols() - 3; i++) // -3避免滑窗出界
      {
        /* evaluate jerk */
        //计算加加速度
        jerk = q.col(i + 3) - 3 * q.col(i + 2) + 3 * q.col(i + 1) - q.col(i); //
        cost += jerk.squaredNorm(); // 加加速度的平方和作为代价
        temp_j = 2.0 * jerk;
        /* jerk gradient */
        //更新梯度
        gradient.col(i + 0) += -temp_j;
        gradient.col(i + 1) += 3.0 * temp_j;
        gradient.col(i + 2) += -3.0 * temp_j;
        gradient.col(i + 3) += temp_j;
      }
    }
    else
    {
      Eigen::Vector3d acc, temp_acc;

      for (int i = 0; i < q.cols() - 2; i++) // -2 避免滑窗出界
      {
        /* evaluate acc */
        acc = q.col(i + 2) - 2 * q.col(i + 1) + q.col(i);
        cost += acc.squaredNorm();
        temp_acc = 2.0 * acc;
        /* acc gradient */
        //更新梯度
        gradient.col(i + 0) += temp_acc;
        gradient.col(i + 1) += -2.0 * temp_acc;
        gradient.col(i + 2) += temp_acc;
      }
    }
  }

  void BsplineOptimizer::calcFeasibilityCost(const Eigen::MatrixXd &q, double &cost,
                                             Eigen::MatrixXd &gradient)
  {

    //#define SECOND_DERIVATIVE_CONTINOUS

#ifdef SECOND_DERIVATIVE_CONTINOUS

    cost = 0.0;
    double demarcation = 1.0; // 1m/s, 1m/s/s
    double ar = 3 * demarcation, br = -3 * pow(demarcation, 2), cr = pow(demarcation, 3);
    double al = ar, bl = -br, cl = cr;

    /* abbreviation */
    double ts, ts_inv2, ts_inv3;
    ts = bspline_interval_;
    ts_inv2 = 1 / ts / ts;
    ts_inv3 = 1 / ts / ts / ts;

    /* velocity feasibility */
    for (int i = 0; i < q.cols() - 1; i++)
    {
      Eigen::Vector3d vi = (q.col(i + 1) - q.col(i)) / ts;

      for (int j = 0; j < 3; j++)
      {
        if (vi(j) > max_vel_ + demarcation)
        {
          double diff = vi(j) - max_vel_;
          cost += (ar * diff * diff + br * diff + cr) * ts_inv3; // multiply ts_inv3 to make vel and acc has similar magnitude

          double grad = (2.0 * ar * diff + br) / ts * ts_inv3;
          gradient(j, i + 0) += -grad;
          gradient(j, i + 1) += grad;
        }
        else if (vi(j) > max_vel_)
        {
          double diff = vi(j) - max_vel_;
          cost += pow(diff, 3) * ts_inv3;
          ;

          double grad = 3 * diff * diff / ts * ts_inv3;
          ;
          gradient(j, i + 0) += -grad;
          gradient(j, i + 1) += grad;
        }
        else if (vi(j) < -(max_vel_ + demarcation))
        {
          double diff = vi(j) + max_vel_;
          cost += (al * diff * diff + bl * diff + cl) * ts_inv3;

          double grad = (2.0 * al * diff + bl) / ts * ts_inv3;
          gradient(j, i + 0) += -grad;
          gradient(j, i + 1) += grad;
        }
        else if (vi(j) < -max_vel_)
        {
          double diff = vi(j) + max_vel_;
          cost += -pow(diff, 3) * ts_inv3;

          double grad = -3 * diff * diff / ts * ts_inv3;
          gradient(j, i + 0) += -grad;
          gradient(j, i + 1) += grad;
        }
        else
        {
          /* nothing happened */
        }
      }
    }

    /* acceleration feasibility */
    for (int i = 0; i < q.cols() - 2; i++)
    {
      Eigen::Vector3d ai = (q.col(i + 2) - 2 * q.col(i + 1) + q.col(i)) * ts_inv2;

      for (int j = 0; j < 3; j++)
      {
        if (ai(j) > max_acc_ + demarcation)
        {
          double diff = ai(j) - max_acc_;
          cost += ar * diff * diff + br * diff + cr;

          double grad = (2.0 * ar * diff + br) * ts_inv2;
          gradient(j, i + 0) += grad;
          gradient(j, i + 1) += -2 * grad;
          gradient(j, i + 2) += grad;
        }
        else if (ai(j) > max_acc_)
        {
          double diff = ai(j) - max_acc_;
          cost += pow(diff, 3);

          double grad = 3 * diff * diff * ts_inv2;
          gradient(j, i + 0) += grad;
          gradient(j, i + 1) += -2 * grad;
          gradient(j, i + 2) += grad;
        }
        else if (ai(j) < -(max_acc_ + demarcation))
        {
          double diff = ai(j) + max_acc_;
          cost += al * diff * diff + bl * diff + cl;

          double grad = (2.0 * al * diff + bl) * ts_inv2;
          gradient(j, i + 0) += grad;
          gradient(j, i + 1) += -2 * grad;
          gradient(j, i + 2) += grad;
        }
        else if (ai(j) < -max_acc_)
        {
          double diff = ai(j) + max_acc_;
          cost += -pow(diff, 3);

          double grad = -3 * diff * diff * ts_inv2;
          gradient(j, i + 0) += grad;
          gradient(j, i + 1) += -2 * grad;
          gradient(j, i + 2) += grad;
        }
        else
        {
          /* nothing happened */
        }
      }
    }

#else

    cost = 0.0;
    /* abbreviation */
    double ts, /*vm2, am2, */ ts_inv2;
    // vm2 = max_vel_ * max_vel_;
    // am2 = max_acc_ * max_acc_;

    ts = bspline_interval_;
    ts_inv2 = 1 / ts / ts;

    /* velocity feasibility */
    for (int i = 0; i < q.cols() - 1; i++)
    {
      Eigen::Vector3d vi = (q.col(i + 1) - q.col(i)) / ts;

      //cout << "temp_v * vi=" ;
      for (int j = 0; j < 3; j++)
      {
        if (vi(j) > max_vel_)
        {
          // cout << "fuck VEL" << endl;
          // cout << vi(j) << endl;
          cost += pow(vi(j) - max_vel_, 2) * ts_inv2; // multiply ts_inv3 to make vel and acc has similar magnitude

          gradient(j, i + 0) += -2 * (vi(j) - max_vel_) / ts * ts_inv2;
          gradient(j, i + 1) += 2 * (vi(j) - max_vel_) / ts * ts_inv2;
        }
        else if (vi(j) < -max_vel_)
        {
          cost += pow(vi(j) + max_vel_, 2) * ts_inv2;

          gradient(j, i + 0) += -2 * (vi(j) + max_vel_) / ts * ts_inv2;
          gradient(j, i + 1) += 2 * (vi(j) + max_vel_) / ts * ts_inv2;
        }
        else
        {
          /* code */
        }
      }
    }

    /* acceleration feasibility */
    for (int i = 0; i < q.cols() - 2; i++)
    {
      Eigen::Vector3d ai = (q.col(i + 2) - 2 * q.col(i + 1) + q.col(i)) * ts_inv2;

      //cout << "temp_a * ai=" ;
      for (int j = 0; j < 3; j++)
      {
        if (ai(j) > max_acc_)
        {
          // cout << "fuck ACC" << endl;
          // cout << ai(j) << endl;
          cost += pow(ai(j) - max_acc_, 2);

          gradient(j, i + 0) += 2 * (ai(j) - max_acc_) * ts_inv2;
          gradient(j, i + 1) += -4 * (ai(j) - max_acc_) * ts_inv2;
          gradient(j, i + 2) += 2 * (ai(j) - max_acc_) * ts_inv2;
        }
        else if (ai(j) < -max_acc_)
        {
          cost += pow(ai(j) + max_acc_, 2);

          gradient(j, i + 0) += 2 * (ai(j) + max_acc_) * ts_inv2;
          gradient(j, i + 1) += -4 * (ai(j) + max_acc_) * ts_inv2;
          gradient(j, i + 2) += 2 * (ai(j) + max_acc_) * ts_inv2;
        }
        else
        {
          /* code */
        }
      }
      //cout << endl;
    }

#endif
  }

  bool BsplineOptimizer::check_collision_and_rebound(void)
  {

    int end_idx = cps_.size - order_;

    /*** Check and segment the initial trajectory according to obstacles ***/
    int in_id, out_id;
    vector<std::pair<int, int>> segment_ids;
    bool flag_new_obs_valid = false;
    int i_end = end_idx - (end_idx - order_) / 3;
    for (int i = order_ - 1; i <= i_end; ++i)
    {

      bool occ = grid_map_->getInflateOccupancy(cps_.points.col(i));

      /*** check if the new collision will be valid ***/
      if (occ)
      {
        for (size_t k = 0; k < cps_.direction[i].size(); ++k)
        {
          cout.precision(2);
          if ((cps_.points.col(i) - cps_.base_point[i][k]).dot(cps_.direction[i][k]) < 1 * grid_map_->getResolution()) // current point is outside all the collision_points.
          {
            occ = false; // Not really takes effect, just for better hunman understanding.
            break;
          }
        }
      }

      if (occ)
      {
        flag_new_obs_valid = true;

        int j;
        for (j = i - 1; j >= 0; --j)
        {
          occ = grid_map_->getInflateOccupancy(cps_.points.col(j));
          if (!occ)
          {
            in_id = j;
            break;
          }
        }
        if (j < 0) // fail to get the obs free point
        {
          printf("ERROR! the drone is in obstacle. This should not happen.");
          in_id = 0;
        }

        for (j = i + 1; j < cps_.size; ++j)
        {
          occ = grid_map_->getInflateOccupancy(cps_.points.col(j));

          if (!occ)
          {
            out_id = j;
            break;
          }
        }
        if (j >= cps_.size) // fail to get the obs free point
        {
          printf("WARN! terminal point of the current trajectory is in obstacle, skip this planning.");

          force_stop_type_ = STOP_FOR_ERROR;
          return false;
        }

        i = j + 1;

        segment_ids.push_back(std::pair<int, int>(in_id, out_id));
      }
    }

    if (flag_new_obs_valid)
    {
      vector<vector<Eigen::Vector3d>> a_star_pathes;
      for (size_t i = 0; i < segment_ids.size(); ++i)
      {
        /*** a star search ***/
        Eigen::Vector3d in(cps_.points.col(segment_ids[i].first)), out(cps_.points.col(segment_ids[i].second));
        if (a_star_->AstarSearch(/*(in-out).norm()/10+0.05*/ 0.1, in, out))
        {
          a_star_pathes.push_back(a_star_->getPath());
        }
        else
        {
          printf("a star error \n");
          segment_ids.erase(segment_ids.begin() + i);
          i--;
        }
      }

      /*** Assign parameters to each segment ***/
      for (size_t i = 0; i < segment_ids.size(); ++i)
      {
        // step 1
        for (int j = segment_ids[i].first; j <= segment_ids[i].second; ++j)
          cps_.flag_temp[j] = false;

        // step 2
        int got_intersection_id = -1;
        for (int j = segment_ids[i].first + 1; j < segment_ids[i].second; ++j)
        {
          Eigen::Vector3d ctrl_pts_law(cps_.points.col(j + 1) - cps_.points.col(j - 1)), intersection_point;
          int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id; // Let "Astar_id = id_of_the_most_far_away_Astar_point" will be better, but it needs more computation
          double val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot(ctrl_pts_law), last_val = val;
          while (Astar_id >= 0 && Astar_id < (int)a_star_pathes[i].size())
          {
            last_Astar_id = Astar_id;

            if (val >= 0)
              --Astar_id;
            else
              ++Astar_id;

            val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot(ctrl_pts_law);

            // cout << val << endl;

            if (val * last_val <= 0 && (abs(val) > 0 || abs(last_val) > 0)) // val = last_val = 0.0 is not allowed
            {
              intersection_point =
                  a_star_pathes[i][Astar_id] +
                  ((a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id]) *
                   (ctrl_pts_law.dot(cps_.points.col(j) - a_star_pathes[i][Astar_id]) / ctrl_pts_law.dot(a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id])) // = t
                  );

              got_intersection_id = j;
              break;
            }
          }

          if (got_intersection_id >= 0)
          {
            cps_.flag_temp[j] = true;
            double length = (intersection_point - cps_.points.col(j)).norm();
            if (length > 1e-5)
            {
              for (double a = length; a >= 0.0; a -= grid_map_->getResolution())
              {
                bool occ = grid_map_->getInflateOccupancy((a / length) * intersection_point + (1 - a / length) * cps_.points.col(j));

                if (occ || a < grid_map_->getResolution())
                {
                  if (occ)
                    a += grid_map_->getResolution();
                  cps_.base_point[j].push_back((a / length) * intersection_point + (1 - a / length) * cps_.points.col(j));
                  cps_.direction[j].push_back((intersection_point - cps_.points.col(j)).normalized());
                  break;
                }
              }
            }
            else
            {
              got_intersection_id = -1;
            }
          }
        }

        //step 3
        if (got_intersection_id >= 0)
        {
          for (int j = got_intersection_id + 1; j <= segment_ids[i].second; ++j)
            if (!cps_.flag_temp[j])
            {
              cps_.base_point[j].push_back(cps_.base_point[j - 1].back());
              cps_.direction[j].push_back(cps_.direction[j - 1].back());
            }

          for (int j = got_intersection_id - 1; j >= segment_ids[i].first; --j)
            if (!cps_.flag_temp[j])
            {
              cps_.base_point[j].push_back(cps_.base_point[j + 1].back());
              cps_.direction[j].push_back(cps_.direction[j + 1].back());
            }
        }
        else
          printf("Failed to generate direction. It doesn't matter.\n");
      }

      force_stop_type_ = STOP_FOR_REBOUND;
      return true;
    }

    return false;
  }

  bool BsplineOptimizer::BsplineOptimizeTrajRebound(Eigen::MatrixXd &optimal_points, double ts)
  {
    setBsplineInterval(ts);

    bool flag_success = rebound_optimize();

    optimal_points = cps_.points;

    return flag_success;
  }

  bool BsplineOptimizer::BsplineOptimizeTrajRefine(const Eigen::MatrixXd &init_points, const double ts, Eigen::MatrixXd &optimal_points)
  {

    setControlPoints(init_points);
    setBsplineInterval(ts);

    bool flag_success = refine_optimize();

    optimal_points = cps_.points;

    return flag_success;
  }

  bool BsplineOptimizer::rebound_optimize()
  {
    iter_num_ = 0;
    int start_id = order_;
    int end_id = this->cps_.size - order_;
    variable_num_ = 3 * (end_id - start_id);
    double final_cost;

   // ros::Time t0 = ros::Time::now(), t1, t2;
    int restart_nums = 0, rebound_times = 0;
    ;
    bool flag_force_return, flag_occ, success;
    new_lambda2_ = lambda2_;
    constexpr int MAX_RESART_NUMS_SET = 3;
    do
    {
      /* ---------- prepare ---------- */
      min_cost_ = std::numeric_limits<double>::max();
      iter_num_ = 0;
      flag_force_return = false;
      flag_occ = false;
      success = false;

      double q[variable_num_];
      memcpy(q, cps_.points.data() + 3 * start_id, variable_num_ * sizeof(q[0]));

      lbfgs::lbfgs_parameter_t lbfgs_params;
      lbfgs::lbfgs_load_default_parameters(&lbfgs_params);
      lbfgs_params.mem_size = 16;
      lbfgs_params.max_iterations = 200;
      lbfgs_params.g_epsilon = 0.01;

      /* ---------- optimize ---------- */
     // t1 = ros::Time::now();
      int result = lbfgs::lbfgs_optimize(variable_num_, q, &final_cost, BsplineOptimizer::costFunctionRebound, NULL, BsplineOptimizer::earlyExit, this, &lbfgs_params);
     // t2 = ros::Time::now();
     // double time_ms = (t2 - t1).toSec() * 1000;
     // double total_time_ms = (t2 - t0).toSec() * 1000;

      /* ---------- success temporary, check collision again ---------- */
      if (result == lbfgs::LBFGS_CONVERGENCE ||
          result == lbfgs::LBFGSERR_MAXIMUMITERATION ||
          result == lbfgs::LBFGS_ALREADY_MINIMIZED ||
          result == lbfgs::LBFGS_STOP)
      {
        //ROS_WARN("Solver error in planning!, return = %s", lbfgs::lbfgs_strerror(result));
        flag_force_return = false;

        UniformBspline traj = UniformBspline(cps_.points, 3, bspline_interval_);
        double tm, tmp;
        traj.getTimeSpan(tm, tmp);
        double t_step = (tmp - tm) / ((traj.evaluateDeBoorT(tmp) - traj.evaluateDeBoorT(tm)).norm() / grid_map_->getResolution());
        for (double t = tm; t < tmp * 2 / 3; t += t_step) // Only check the closest 2/3 partition of the whole trajectory.
        {
          flag_occ = grid_map_->getInflateOccupancy(traj.evaluateDeBoorT(t));
          if (flag_occ)
          {
            //cout << "hit_obs, t=" << t << " P=" << traj.evaluateDeBoorT(t).transpose() << endl;

            if (t <= bspline_interval_) // First 3 control points in obstacles!
            {
              cout << cps_.points.col(1).transpose() << "\n"
                   << cps_.points.col(2).transpose() << "\n"
                   << cps_.points.col(3).transpose() << "\n"
                   << cps_.points.col(4).transpose() << endl;
              printf("First 3 control points in obstacles! return false, t=%f", t);
              return false;
            }

            break;
          }
        }

        if (!flag_occ)
        {
          //printf("\033[32miter(+1)=%d,time(ms)=%5.3f,total_t(ms)=%5.3f,cost=%5.3f\n\033[0m", iter_num_, time_ms, total_time_ms, final_cost);
          success = true;
        }
        else // restart
        {
          restart_nums++;
          initControlPoints(cps_.points, false);
          new_lambda2_ *= 2;

          //printf("\033[32miter(+1)=%d,time(ms)=%5.3f,keep optimizing\n\033[0m", iter_num_, time_ms);
        }
      }
      else if (result == lbfgs::LBFGSERR_CANCELED)
      {
        flag_force_return = true;
        rebound_times++;
        cout << "iter=" << iter_num_ << ",rebound." << endl;
      }
      else
      {
        printf("Solver error. Return = %d, %s. Skip this planning.\n", result, lbfgs::lbfgs_strerror(result));
        // while (ros::ok());
      }

    } while ((flag_occ && restart_nums < MAX_RESART_NUMS_SET) ||
             (flag_force_return && force_stop_type_ == STOP_FOR_REBOUND && rebound_times <= 20));

    return success;
  }

  //执行优化
  bool BsplineOptimizer::refine_optimize()
  {
    //记录迭代次数
    iter_num_ = 0;
    //设置范围（起始点，终止点）
    int start_id = order_;
    int end_id = this->cps_.points.cols() - order_;
    //变量数
    //3 * 点的个数（每个点3个double
    variable_num_ = 3 * (end_id - start_id);

    // 申请空间
    double q[variable_num_];
    double final_cost; //记录最终代价

    //复制数据到q中
    memcpy(q, cps_.points.data() + 3 * start_id, variable_num_ * sizeof(q[0]));

    double origin_lambda4 = lambda4_;
    bool flag_safe = true;
    int iter_count = 0;
    do
    {
      //初始化lbfgs参数
      lbfgs::lbfgs_parameter_t lbfgs_params;
      lbfgs::lbfgs_load_default_parameters(&lbfgs_params);
      lbfgs_params.mem_size = 16;
      lbfgs_params.max_iterations = 200;
      lbfgs_params.g_epsilon = 0.001;
      //传入变量和代价函数执行优化
      int result = lbfgs::lbfgs_optimize(variable_num_, q, &final_cost, BsplineOptimizer::costFunctionRefine, NULL, NULL, this, &lbfgs_params);
      //正常结束优化
      if (result == lbfgs::LBFGS_CONVERGENCE ||
          result == lbfgs::LBFGSERR_MAXIMUMITERATION ||
          result == lbfgs::LBFGS_ALREADY_MINIMIZED ||
          result == lbfgs::LBFGS_STOP)
      {
        //pass
      }
      else
      {
        printf("Solver error in refining!, return = %d, %s", result, lbfgs::lbfgs_strerror(result));
      }

      //对优化完的控制点构造曲线检查
      UniformBspline traj = UniformBspline(cps_.points, 3, bspline_interval_);
      //
      double tm, tmp;
      traj.getTimeSpan(tm, tmp);
      double t_step = (tmp - tm) / ((traj.evaluateDeBoorT(tmp) - traj.evaluateDeBoorT(tm)).norm() / grid_map_->getResolution()); // Step size is defined as the maximum size that can passes throgth every gird.
      for (double t = tm; t < tmp * 2 / 3; t += t_step)
      {
        //如果碰到障碍物，则返回失败
        if (grid_map_->getInflateOccupancy(traj.evaluateDeBoorT(t)))
        {
          // cout << "Refined traj hit_obs, t=" << t << " P=" << traj.evaluateDeBoorT(t).transpose() << endl;
          Eigen::MatrixXd ref_pts(ref_pts_.size(), 3);
          for (size_t i = 0; i < ref_pts_.size(); i++)
          {
            ref_pts.row(i) = ref_pts_[i].transpose();
          }

          flag_safe = false;
          break;
        }
      }

      if (!flag_safe)
        lambda4_ *= 2;

      iter_count++;
    } while (!flag_safe && iter_count <= 0);

    lambda4_ = origin_lambda4;

    //cout << "iter_num_=" << iter_num_ << endl;

    return flag_safe;
  }

  // 合并代价函数1: 
  void BsplineOptimizer::combineCostRebound(const double *x, double *grad, double &f_combine, const int n)
  {

    memcpy(cps_.points.data() + 3 * order_, x, n * sizeof(x[0])); //复制数据

    /* ---------- evaluate cost and gradient ---------- */
    double f_smoothness, f_distance, f_feasibility;

    Eigen::MatrixXd g_smoothness = Eigen::MatrixXd::Zero(3, cps_.size);
    Eigen::MatrixXd g_distance = Eigen::MatrixXd::Zero(3, cps_.size);
    Eigen::MatrixXd g_feasibility = Eigen::MatrixXd::Zero(3, cps_.size);

    // 平滑代价
    calcSmoothnessCost(cps_.points, f_smoothness, g_smoothness);
    // 障碍物距离代价
    calcDistanceCostRebound(cps_.points, f_distance, g_distance, iter_num_, f_smoothness);
    // 运动学可行性代价
    calcFeasibilityCost(cps_.points, f_feasibility, g_feasibility);

    f_combine = lambda1_ * f_smoothness + new_lambda2_ * f_distance + lambda3_ * f_feasibility;
    //printf("origin %f %f %f %f\n", f_smoothness, f_distance, f_feasibility, f_combine);

    Eigen::MatrixXd grad_3D = lambda1_ * g_smoothness + new_lambda2_ * g_distance + lambda3_ * g_feasibility;
    memcpy(grad, grad_3D.data() + 3 * order_, n * sizeof(grad[0]));
  }

  // 合并代价函数2: 
  void BsplineOptimizer::combineCostRefine(const double *x, double *grad, double &f_combine, const int n)
  {

    memcpy(cps_.points.data() + 3 * order_, x, n * sizeof(x[0]));

    /* ---------- evaluate cost and gradient ---------- */
    double f_smoothness, f_fitness, f_feasibility;

    Eigen::MatrixXd g_smoothness = Eigen::MatrixXd::Zero(3, cps_.points.cols());
    Eigen::MatrixXd g_fitness = Eigen::MatrixXd::Zero(3, cps_.points.cols());
    Eigen::MatrixXd g_feasibility = Eigen::MatrixXd::Zero(3, cps_.points.cols());

    //time_satrt = ros::Time::now();

    // 平滑代价
    calcSmoothnessCost(cps_.points, f_smoothness, g_smoothness);
    // 
    calcFitnessCost(cps_.points, f_fitness, g_fitness);
    // 运动学可行性代价
    calcFeasibilityCost(cps_.points, f_feasibility, g_feasibility);

    /* ---------- convert to solver format...---------- */
    f_combine = lambda1_ * f_smoothness + lambda4_ * f_fitness + lambda3_ * f_feasibility;
    // printf("origin %f %f %f %f\n", f_smoothness, f_fitness, f_feasibility, f_combine);

    Eigen::MatrixXd grad_3D = lambda1_ * g_smoothness + lambda4_ * g_fitness + lambda3_ * g_feasibility;
    memcpy(grad, grad_3D.data() + 3 * order_, n * sizeof(grad[0]));
  }

} // namespace ego_planner