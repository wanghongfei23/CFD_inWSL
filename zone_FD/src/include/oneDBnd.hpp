#pragma once

#include "data.hpp"
#include "info.hpp"
#include "macro.hpp"
class OneDBnd {
public:
  OneDBnd() {};
  OneDBnd(int, int, BndType);
  real &operator()(int, int);
  std::vector<real>::iterator operator()(int);

  int getN();
  BndType getType();
  void setUpdate(Data *, int, int);
  void setValue(std::vector<real>);
  void setInfo(Info *);
  void setCoor(std::array<real, 3>, std::array<real, 3>);

  void update();

private:
  BndType type = TYPENULL;
  std::vector<real> data;
  Data *prim;
  int i0 = 0, offset = 1;

  int n = 0;
  int nVar = 0;

  // for Double Mach problem
  Info *info = nullptr;
  std::array<real, 3> coor = {0, 0, 0}, dh = {0, 0, 0};
};