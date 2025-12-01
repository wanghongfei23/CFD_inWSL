#pragma once
#include "block.hpp"
#include "cgnslib.h"
#include "data.hpp"
#include "info.hpp"
#include "macro.hpp"

// CGNS输入输出类，用于处理CGNS格式的网格和解文件
class CgnsIO {
public:
  // 网格和解数据输出函数
  void BlockCgnsOutput(Block *block, Info *info);
  void solCgnsOutput(Data *data, Info *info);
  void oneDsolOutput(Info *info);

private:
};
