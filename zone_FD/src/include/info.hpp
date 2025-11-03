#pragma once
#include "macro.hpp"
// Info类用于存储和管理计算流体力学(CFD)模拟中的各种参数和配置信息
class Info {
public:
  // 枚举类型 成员变量 = 枚举值

  // 基本求解器参数
  InterMethod spMethod = WCNS5;                         // 插值方法
  EquationType eqType = EULER;                          // 控制方程类型
  int nCase = 2;                                        // 案例编号
  real t = 0;                                           // 当前时间
  real dt = 0.0001;                                     // 时间步长
  int step = 0;                                         // 当前步数
  int endStep = 25000;                                  // 结束步数
  int dim = 1;                                          // 空间维度

  // 输出控制参数
  int outputInterval = 100;                             // 输出间隔步数
  real outputDt = 0.01;                                 // 输出时间间隔
  real outputT = 0;                                     // 上次输出时间

  // // 添加重力加速度参数（RT不稳定性需要）
  // real gravity = 0.0; 

  // 时间步长控制参数
  real CFL = 0.5;                                       // CFL数
  bool fixedtimeSteps = true;                           // 是否使用固定时间步长

  // 隐式求解器参数
  real implicitCFL = 0.01;                              // 隐式求解器的CFL数
  real maxImplicitStep = 100;                           // 最大隐式步数

  // 数值方法参数
  DiffMethod diffMethod = MND6;                         // 差分方法
  InterMethod interMethod = WCNS5;                      // 插值方法，只影响插值权
  SourceType sourceType = SOURCENULL;                   // 源项类型
  FluxMethod fluxMethod = ROE;                          // 通量计算方法

  // 网格参数
  std::array<int, 3> iMax{201, 201, 2};                // 网格最大索引
  std::array<double, 6> calZone{0, 0.3, 0, 0.3, 0, 2}; // 计算区域范围

  // 构造函数和成员函数
  Info();                                               // 构造函数
  int nGhostCell();                                     // 获取Ghost Cell数量
  int nFluxPoint();                                     // 获取通量点数量
  int nPrim();                                          // 获取基本变量数量
  int nCons();                                          // 获取守恒变量数量
  int getDim();                                         // 获取计算域维度

  // 网格相关参数和函数
  bool constH = true;                                   // 是否使用恒定网格步长
  real interval = 0;                                    // 网格间隔
  real geth(int);                                       // 获取指定方向上的网格步长

  // 边界条件和文件操作相关函数
  BndType defaultBndType();                             // 获取默认边界类型
  std::array<int, 3> icMax();                           // 获取内部网格点最大索引
  std::string filename();                               // 获取输出文件名

  // 变量名称列表获取函数
  std::vector<std::string> getVarNameListCons();        // 获取守恒变量名称列表
  std::vector<std::string> getVarNameListPrim();        // 获取基本变量名称列表
  std::vector<std::string> getVarNameListRhs();         // 获取RHS变量名称列表
};
