#include "oneDBnd.hpp"

OneDBnd::OneDBnd(int n_, int nVar_, BndType bType_) {
  n = n_;
  nVar = nVar_;
  data.resize(n * nVar, 0.0);
  type = bType_;
}
real &OneDBnd::operator()(int i, int ivar) { return data.at(i * nVar + ivar); }

std::vector<real>::iterator OneDBnd::operator()(int i) {
  return data.begin() + i * nVar;
}

int OneDBnd::getN() { return n; }

void OneDBnd::setValue(std::vector<real> value) {
  if (value.size() != data.size()) {
    std::cout << "OneDBnd error: setValue() incorrect size\n";
  }
  std::copy(value.begin(), value.end(), data.begin());
}
BndType OneDBnd::getType() { return type; }

void OneDBnd::setUpdate(Data *prim_, int i0_, int offset_) {
  prim = prim_;
  i0 = i0_;
  offset = offset_;
}

void OneDBnd::update() {
  switch (type) {
  case PERIODIC1D:
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nVar; j++) {
        data[i * nVar + j] = (*prim)(i0 + i * offset, j);
      }
    }
    break;
  case SUPERSONICOUTLET:
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nVar; j++) {
        data[i * nVar + j] = (*prim)(i0, j);
      }
    }
    break;
  case SYMMETRYX:
    // only for 2D
    for (int i = 0; i < n; i++) {
      // data[i*nVar+0]=(*prim)(i0,0);
      // data[i*nVar+1]=-(*prim)(i0,1);
      // data[i*nVar+2]=(*prim)(i0,2);
      // data[i*nVar+3]=(*prim)(i0,3);

      data[i * nVar + 0] = (*prim)(i0 + i * offset, 0);
      data[i * nVar + 1] = -(*prim)(i0 + i * offset, 1);
      data[i * nVar + 2] = (*prim)(i0 + i * offset, 2);
      data[i * nVar + 3] = (*prim)(i0 + i * offset, 3);
    }
    break;
  case SYMMETRY1D:
    for (int i = 0; i < n; i++) {
      // data[i*nVar+0]=(*prim)(i0,0);
      // data[i*nVar+1]=-(*prim)(i0,1);
      // data[i*nVar+2]=(*prim)(i0,2);
      // data[i*nVar+3]=(*prim)(i0,3);

      data[i * nVar + 0] = (*prim)(i0 + i * offset, 0);
      data[i * nVar + 1] = -(*prim)(i0 + i * offset, 1);
      data[i * nVar + 2] = (*prim)(i0 + i * offset, 2);
    }
    break;
  case SYMMETRYY:
    // only for 2D
    for (int i = 0; i < n; i++) {
      data[i * nVar + 0] = (*prim)(i0 + i * offset, 0);
      data[i * nVar + 1] = (*prim)(i0 + i * offset, 1);
      data[i * nVar + 2] = -(*prim)(i0 + i * offset, 2);
      data[i * nVar + 3] = (*prim)(i0 + i * offset, 3);

      // data[i*nVar+0]=(*prim)(i0,0);
      // data[i*nVar+1]=(*prim)(i0,1);
      // data[i*nVar+2]=-(*prim)(i0,2);
      // data[i*nVar+3]=(*prim)(i0,3);
    }
    break;
  case DoubleMachUp:
    // only for 2D
    {
      for (int i = 0; i < n; i++) {
        real x = coor[0];
        real y = -dh[1] / 2 + i * dh[1];
        // real y=dh[1]/2;
        real gt = 1.0 / 6.0 + std::sqrt(3.0) / 3.0 * (1.0 + 20 * info->t);
        std::array<real, 4> exactValues;
        if (y > std::sqrt(3.0) * (x - gt)) {
          exactValues = {8.0, 8.25 * cos(M_PI / 6), -8.25 * sin(M_PI / 6),
                         116.5};
        } else {
          exactValues = {1.4, 0, 0, 1.0};
        }

        data[i * nVar + 0] = exactValues[0];
        data[i * nVar + 1] = exactValues[1];
        data[i * nVar + 2] = exactValues[2];
        data[i * nVar + 3] = exactValues[3];
      }
    }
    break;

  default:
    break;
  }
}

void OneDBnd::setInfo(Info *info_) { info = info_; }
void OneDBnd::setCoor(std::array<real, 3> coor_, std::array<real, 3> dh_) {
  coor = coor_;
  dh = dh_;
}