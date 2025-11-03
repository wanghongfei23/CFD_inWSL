#pragma once
#include "data.hpp"
#include "info.hpp"
#include "oneDBnd.hpp"

// 边界条件类，用于管理所有边界条件
class Bnds {
public:
  // 边界条件获取和更新函数
  std::array<std::shared_ptr<OneDBnd>, 2> getOneDBnd(int, int, int);
  void update();

private:
  // 允许初始化器访问私有成员
  friend class Initializer;
  
  // 边界条件相关数据
  int dim;                                                    // 问题维度
  std::array<int, 3> iMax;                                    // 网格最大索引
  std::vector<std::shared_ptr<OneDBnd>> oneDBnds;             // 一维边界条件数组，包含所有方向的边界条件
};