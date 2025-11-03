#include "eigenSystem.hpp"

eigensystemEuler2D::eigensystemEuler2D(const std::array<real, 4> &prim,
                                       const std::array<real, 3> &norm_) {
  r = prim[0], u = prim[1], v = prim[2], p = prim[3];
  gamma = GAMMA, ek = (u * u + v * v) / 2;
  h = p / r * gamma / (1 - gamma);
  c = std::sqrt(gamma * p / r);
  norm = norm_;
  Vn = norm[0] * u + norm[1] * v;
}

eigensystemEuler2D::eigensystemEuler2D(const std::array<real, 4> &priml,
                                       const std::array<real, 4> &primr,
                                       const std::array<real, 3> &norm_) {
  //chenyuqing: 构建左右特征变换矩阵
  norm = norm_;
 
  gamma = GAMMA;
  enum { LL, RR };
  real rl = priml[0], ul = priml[1], vl = priml[2], pl = priml[3];
  real rr = primr[0], ur = primr[1], vr = primr[2], pr = primr[3];

  std::array<real, 2> H;
  H[LL] = (ul * ul + vl * vl) / 2 + pl / rl * gamma / (gamma - 1);
  H[RR] = (ur * ur + vr * vr) / 2 + pr / rr * gamma / (gamma - 1);
  real coef1 = std::sqrt(rl);
  real coef2 = std::sqrt(rr);
  real divisor = 1.0 / (coef1 + coef2);

  r = std::sqrt(rl * rr);
  u = (coef1 * ul + coef2 * ur) * divisor;
  v = (coef1 * vl + coef2 * vr) * divisor;
  real ht = (coef1 * H[LL] + coef2 * H[RR]) * divisor;
  Vn = norm[0] * u + norm[1] * v;
  ek = (u * u + v * v) / 2;
  h = ht -ek;
  //chenyuqing：修改特征矩阵Roe平均值
  // real co1 = std::sqrt(rl) / (std::sqrt(rl) + std::sqrt(rr));
  // real co2 = 1.0 - co1;
  // r = std::sqrt(rl * rr);
  // u = co1 * ul + co2 * ur;
  // v = co1 * vl + co2 * vr;
  // real ht = co1 * H[LL] + co2 * H[RR];
  // Vn = norm[0] * u + norm[1] * v;
  // ek = (u * u + v * v) / 2;
  // h = ht - u * u / 2 - v * v / 2;


  c = std::sqrt((gamma - 1) * h);
  p = r * h * ((gamma - 1) / gamma);

  leftEig = {(-norm[Y] * u + norm[X] * v),
             norm[Y],
             (-norm[X]),
             0,
             (h - ek),
             u,
             v,
             -1.0,
             (Vn / c + ek / h) / 2,
             (-norm[X] / c - u / h) / 2,
             (-norm[Y] / c - v / h) / 2,
             1.0 / (2.0 * h),
             (-Vn / c + ek / h) / 2,
             (norm[X] / c - u / h) / 2,
             (norm[Y] / c - v / h) / 2,
             1.0 / (2.0 * h)};

  rightEig = {0,
              1.0 / h,
              1,
              1,
              norm[Y],
              u / h,
              u - norm[X] * c,
              u + norm[X] * c,
              -norm[X],
              v / h,
              v - norm[Y] * c,
              v + norm[Y] * c,
              norm[Y] * u - norm[X] * v,
              ek / h,
              h + ek - Vn * c,
              h + ek + Vn * c};
}

std::array<real, 4>
eigensystemEuler2D::primToChar(const std::array<real, 4> &prim) {
  //chenyuqing: 原始变量转为特征变量
  real rt = prim[0], ut = prim[1], vt = prim[2], pt = prim[3];
  real ekt = (ut * ut + vt * vt) / 2;
  real rut = rt * ut, rvt = rt * vt, ret = pt / (gamma - 1) + rt * ekt;

  std::array<real, 4> res;

  res[0] =
      rut * leftEig[1] + rvt * leftEig[2] + ret * leftEig[3] + rt * leftEig[0];

  res[1] =
      rut * leftEig[5] + rvt * leftEig[6] + ret * leftEig[7] + rt * leftEig[4];

  res[2] = rut * leftEig[9] + rvt * leftEig[10] + ret * leftEig[11] +
           rt * leftEig[8];

  res[3] = rut * leftEig[13] + rvt * leftEig[14] + ret * leftEig[15] +
           rt * leftEig[12];

  return res;
  //chenyuqing：修改原始变量转换成特征变量
  // res[0] =
  //     rt * leftEig[0] + rut * leftEig[1] + rvt * leftEig[2] + ret * leftEig[3];

  // res[1] =
  //     rt * leftEig[4] + rut * leftEig[5] + rvt * leftEig[6] + ret * leftEig[7] ;

  // res[2] = rt * leftEig[8] + rut * leftEig[9] + rvt * leftEig[10] + ret * leftEig[11] ;

  // res[3] = rt * leftEig[12] + rut * leftEig[13] + rvt * leftEig[14] + ret * leftEig[15] ;
       
  // return res;

}

std::array<real, 4>
eigensystemEuler2D::primToCons(const std::array<real, 4> &prim) {
  //chenyuqing: 原始变量转为守恒变量
  real rt = prim[0], ut = prim[1], vt = prim[2], pt = prim[3];
  real ekt = (ut * ut + vt * vt) / 2;
  real rut = rt * ut, rvt = rt * vt, ret = pt / (gamma - 1) + rt * ekt;

  return {rt,rut,rvt,ret};
}
std::array<real, 4>
eigensystemEuler2D::charToPrim(const std::array<real, 4> &chars) {
  //chenyuqing: 特征变量转为原始变量
  real ch1 = chars[0], ch2 = chars[1], ch3 = chars[2], ch4 = chars[3], rt, rut,
       rvt, ret;
  rt = ch3 * rightEig[2] + ch4 * rightEig[3] + ch2 * rightEig[1] +
       ch1 * rightEig[0];

  rut = ch3 * rightEig[6] + ch4 * rightEig[7] + ch2 * rightEig[5] +
        ch1 * rightEig[4];

  rvt = ch3 * rightEig[10] + ch4 * rightEig[11] + ch2 * rightEig[9] +
        ch1 * rightEig[8];

  ret = ch3 * rightEig[14] + ch4 * rightEig[15] + ch2 * rightEig[13] +
        ch1 * rightEig[12];
  //chenyuqing：修改特征变量转为原始变量
  // real ch1 = chars[0], ch2 = chars[1], ch3 = chars[2], ch4 = chars[3], rt, rut,
  //      rvt, ret;
  // rt =  ch1 * rightEig[0] + ch2 * rightEig[1] + ch3 * rightEig[2] + ch4 * rightEig[3] ;

  // rut = ch1 * rightEig[4] + ch2 * rightEig[5] + ch3 * rightEig[6] + ch4 * rightEig[7] 
  //       ;

  // rvt =  ch1 * rightEig[8] + ch2 * rightEig[9] + ch3 * rightEig[10] + ch4 * rightEig[11] 
  //      ;

  // ret =  ch1 * rightEig[12] + ch2 * rightEig[13] + ch3 * rightEig[14] + ch4 * rightEig[15] ;

  real ut = rut / rt;
  real vt = rvt / rt;
  real rekt = (rut * rut + rvt * rvt) / rt / 2;
  real pt = (gamma - 1) * (ret - rekt);
  return {rt, ut, vt, pt};
}

std::array<real, 4>
eigensystemEuler2D::charToCons(const std::array<real, 4> &chars) {
  //chenyuqing: 特征变量转为守恒变量
  real ch1 = chars[0], ch2 = chars[1], ch3 = chars[2], ch4 = chars[3], rt, rut,
       rvt, ret;
  rt = ch3 * rightEig[2] + ch4 * rightEig[3] + ch2 * rightEig[1] +
       ch1 * rightEig[0];

  rut = ch3 * rightEig[6] + ch4 * rightEig[7] + ch2 * rightEig[5] +
        ch1 * rightEig[4];

  rvt = ch3 * rightEig[10] + ch4 * rightEig[11] + ch2 * rightEig[9] +
        ch1 * rightEig[8];

  ret = ch3 * rightEig[14] + ch4 * rightEig[15] + ch2 * rightEig[13] +
        ch1 * rightEig[12];
  return {rt, rut, rvt, ret};
  //chenyuqing：修改特征变量转守恒变量
  // real ch1 = chars[0], ch2 = chars[1], ch3 = chars[2], ch4 = chars[3], rt, rut,
  // rvt, ret;
  // rt = ch1 * rightEig[0] +  ch4 * rightEig[3] + ch2 * rightEig[1] + ch3 * rightEig[2] ;

  // rut = ch1 * rightEig[4] + ch4 * rightEig[7] + ch2 * rightEig[5] + ch3 * rightEig[6] ;

  // rvt = ch1 * rightEig[8] + ch4 * rightEig[11] + ch2 * rightEig[9] + ch3 * rightEig[10] ;

  // ret = ch1 * rightEig[12] + ch4 * rightEig[15] + ch2 * rightEig[13] + ch3 * rightEig[14] ;

  // return {rt, rut, rvt, ret};
}

std::array<real, 4>
eigensystemEuler2D::consToPrim(const std::array<real, 4> &cons) {
  //chenyuqing: 守恒变量转为原始变量
  real rt=cons[0],rut=cons[1],rvt=cons[2],ret=cons[3];

    real ut = rut / rt;
    real vt = rvt / rt;
    real rekt = (rut * rut + rvt * rvt) / rt / 2;
    real pt = (gamma - 1) * (ret - rekt);
    //chenyuqing：修改守恒变量转为原始变量
    //real pt = (gamma - 1) * (ret - (rut * rut) / rt / 2 - rvt * rvt / rt / 2);
    return {rt, ut, vt, pt};
}



eigensystemEuler1D::eigensystemEuler1D(const std::array<real, 3> &priml,
                                       const std::array<real, 3> &primr) {
  gamma = GAMMA;
  enum { LL, RR };
  real rl = priml[0], ul = priml[1], pl = priml[2];
  real rr = primr[0], ur = primr[1], pr = primr[2];

  std::array<real, 2> H;
  H[LL] = (ul * ul) / 2 + pl / rl * gamma / (gamma - 1);
  H[RR] = (ur * ur) / 2 + pr / rr * gamma / (gamma - 1);
  real coef1 = std::sqrt(rl);
  real coef2 = std::sqrt(rr);
  real divisor = 1.0 / (std::sqrt(rl) + std::sqrt(rr));

  r = std::sqrt(rl * rr);
  u = (coef1 * ul + coef2 * ur) * divisor;
  real ht = (coef1 * H[LL] + coef2 * H[RR]) * divisor;
  ek = (u * u) / 2;
  h = ht - ek;
  c = std::sqrt((gamma - 1) * h);
  p = r * h * ((gamma - 1) / gamma);

  real gamma_1 = gamma - 1;
  leftEig = {ht + c * (u - c) / gamma_1,        -u - c / gamma_1, 1.0,
             -2.0 * ht + 4.0 * c * c / gamma_1, 2.0 * u,          -2.0,
             ht - c * (u + c) / gamma_1,        -u + c / gamma_1, 1.0};
  real factorEig = 0.5 * gamma_1 / (c * c);

  for (int ii = 0; ii < 9; ii++)
    leftEig[ii] *= factorEig;

  rightEig = {1.0,   1.0,        1.0,         u - c,     u,
              u + c, ht - u * c, 0.5 * u * u, ht + u * c};
}

std::array<real, 3>
eigensystemEuler1D::primToChar(const std::array<real, 3> &prim) {
  real rt = prim[0], ut = prim[1], pt = prim[2];
  real ekt = (ut * ut) / 2;
  real rut = rt * ut, ret = pt / (gamma - 1) + rt * ekt;

  std::array<real, 3> res;

  res[0] = rt * leftEig[0] + rut * leftEig[1] + ret * leftEig[2];

  res[1] = rt * leftEig[3] + rut * leftEig[4] + ret * leftEig[5];

  res[2] = rt * leftEig[6] + rut * leftEig[7] + ret * leftEig[8];

  return res;
}

std::array<real, 3>
eigensystemEuler1D::primToCons(const std::array<real, 3> &prim) {
  real rt = prim[0], ut = prim[1], pt = prim[2];
  real ekt = (ut * ut) / 2;
  real rut = rt * ut, ret = pt / (gamma - 1) + rt * ekt;
  return {rt,rut,ret};
}

std::array<real, 3>
eigensystemEuler1D::charToPrim(const std::array<real, 3> &chars) {
  real ch1 = chars[0], ch2 = chars[1], ch3 = chars[2], rt, rut, ret;

  rt = ch1 * rightEig[0] + ch2 * rightEig[1] + ch3 * rightEig[2];

  rut = ch1 * rightEig[3] + ch2 * rightEig[4] + ch3 * rightEig[5];

  ret = ch1 * rightEig[6] + ch2 * rightEig[7] + ch3 * rightEig[8];

  real ut = rut / rt;
  real Et = ret / rt;
  real ekt = (ut * ut) / 2;
  real pt = (gamma - 1) * (ret - rt * ekt);
  return {rt, ut, pt};
}

std::array<real, 3>
eigensystemEuler1D::charToCons(const std::array<real, 3> &chars) {
  // 特征变量转为守恒变量
  real ch1 = chars[0], ch2 = chars[1], ch3 = chars[2];

  // 计算守恒变量
  real rt = ch1 * rightEig[0] + ch2 * rightEig[1] + ch3 * rightEig[2]; // 质量 (rho)
  real rut = ch1 * rightEig[3] + ch2 * rightEig[4] + ch3 * rightEig[5]; // 动量 (rho * u)
  real ret = ch1 * rightEig[6] + ch2 * rightEig[7] + ch3 * rightEig[8]; // 总能量 (rho * E)

  // 返回守恒变量
  return {rt, rut, ret};
}

std::array<real, 3>
eigensystemEuler1D::consToPrim(const std::array<real, 3> &cons) {
  // 守恒变量转为原始变量
  real rt = cons[0]; // 质量 (rho)
  real rut = cons[1]; // 动量 (rho * u)
  real ret = cons[2]; // 总能量 (rho * E)

  // 计算原始变量
  real ut = rut / rt; // 速度 (u)
  real Et = ret / rt; // 单位质量总能量 (E)
  real ekt = (ut * ut) / 2; // 动能 (u^2 / 2)
  real pt = (gamma - 1) * (ret - rt * ekt); // 压力 (p)

  // 返回原始变量
  return {rt, ut, pt};
}