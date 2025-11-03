#include "macro.hpp"

// 差分格式基类，定义差分格式的接口
class DifferenceScheme {
public:
  // 差分格式所需的数据点数量
  const int nFluxHalf = 0;
  const int nFluxNode = 0;
  
  // 差分计算接口函数
  virtual real operator()(std::vector<real>::iterator fluxHalf,
                          std::vector<real>::iterator fluxNode, int nvar) = 0;
};

// 二阶差分格式类
class DifSecondOrder : public DifferenceScheme {
public:
  // 差分格式所需的数据点数量
  const int nFluxHalf = 0;
  const int nFluxNode = 0;
  
  // 二阶差分计算实现
  real operator()(std::vector<real>::iterator fluxHalf,
                  std::vector<real>::iterator fluxNode, int nvar) override {
    return fluxHalf[1] - fluxHalf[0];
  };
};

// 二阶差分格式函数
real secondOrder(std::vector<real>::iterator fluxHalf,
                 std::vector<real>::iterator fluxNode, int nvar) {
  return fluxHalf[1] - fluxHalf[0];
}

// 显式六阶差分格式函数
real explicit6(std::vector<real>::iterator fluxHalf,
               std::vector<real>::iterator fluxNode, int nvar) {
  constexpr std::array<real, 3> w = {75.0 / 64.0, 25.0 / (128.0 * 3.0),
                                     3.0 / (128.0 * 50)};
  return w[0] * (*(fluxHalf - nvar * 3) - *(fluxHalf - nvar * 2)) -
         w[1] * (*(fluxHalf - nvar * 4) - *(fluxHalf - nvar * 1)) +
         w[2] * (*(fluxHalf - nvar * 5) - *(fluxHalf - nvar * 0));
}

// // 六阶单调性保持差分格式函数
// real MND6(std::vector<real>::iterator fluxHalf,
//           std::vector<real>::iterator fluxNode, int nvar) {
//   constexpr std::array<real, 3> w = {3.0 / 2.0, -3.0 / 10.0, 1.0 / 30.0};
//   return w[1] * (*(fluxNode - nvar * 1) - *(fluxNode - nvar * 0)) +
//          w[0] * (*(fluxHalf - nvar * 2) - *(fluxHalf - nvar * 1)) +
//          w[2] * (*(fluxHalf - nvar * 3) - *(fluxHalf - nvar * 0));
// }