#pragma once
#include "data.hpp"
#include "macro.hpp"
#include "oneDBnd.hpp"

// 数据操作器基类，用于处理数据的边界条件和求解操作
class DataManipulater {
public:
  // 构造函数
  DataManipulater() {};
  DataManipulater(DataReader varsReader_, std::shared_ptr<OneDBnd> bndL_,
                  std::shared_ptr<OneDBnd> bndR_);
  
  // 初始化函数
  void init(DataReader varsReader_, std::shared_ptr<OneDBnd> bndL_,
            std::shared_ptr<OneDBnd> bndR_);
  
  // 设置法向量函数
  void setConstNorm(std::array<real, 3> norm_);
  
  // 数据获取函数
  std::shared_ptr<Data> getData() { return this->data; }
  
  // 求解接口函数
  virtual void solve() = 0;
  
  // 测试函数
  void test() { std::cout << "test success\n"; }

  // 数据指针
  std::shared_ptr<Data> data;

protected:
  // 迭代器
  std::vector<real>::iterator iter;
  
  // 纯虚函数接口
  virtual void initVarsR() = 0;
  virtual void leftBnd() = 0;
  virtual void internal() = 0;
  virtual void rightBnd() = 0;

  // 数据读取和边界条件
  DataReader varsReader;
  std::shared_ptr<OneDBnd> bndL, bndR;
  
  // 法向量和网格参数
  std::array<real, 3> norm;
  int n, nvar;
};