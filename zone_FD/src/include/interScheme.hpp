#pragma once
#include "macro.hpp"
#include <array>
// inline real weno5_JSchen(std::array<real, 5>);
inline real u1(real q1, real q2, real q3) {
  return 3.0 / 8.0 * q1 - 5.0 / 4.0 * q2 + 15.0 / 8.0 * q3;
}
inline real u2(real q1, real q2, real q3) {
  return -1.0 / 8.0 * q1 + 3.0 / 4.0 * q2 + 3.0 / 8.0 * q3;
}
inline real u3(real q1, real q2, real q3) {
  return 3.0 / 8.0 * q1 + 3.0 / 4.0 * q2 - 1.0 / 8.0 * q3;
}

inline real weno5_JSchen(std::array<real, 5> q) {
  real eps = 1e-6;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

  real sumbeta = 0, result = 0;
  for (int i = 0; i < 3; i++) {
    beta[i] = gamma[i] / pow(eps + beta[i], 2);
    sumbeta += beta[i];
  }
  for (int i = 0; i < 3; i++)
    result += beta[i] * u[i];
  return result / sumbeta;
}

inline real weno5_Cong(std::array<real, 5> q) {
  real eps = 1e-14;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta, u;
  // beta[0]= 0.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
  //          + 1.0/1.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

  //  beta[1]= 0.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
  //          + 1.0/1.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

  //  beta[2]= 0.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
  //          + 1.0/1.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

  beta[0] =
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) +
      10 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2) / (q[2] * q[2]);

  beta[1] =
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) +
      10 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2) / (q[2] * q[2]);

  beta[2] =
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) +
      10 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2) / (q[2] * q[2]);

  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  real sumbeta = 0, result = 0;
  for (int i = 0; i < 3; i++) {
    beta[i] = gamma[i] / pow(eps + beta[i], 2.0);
    sumbeta += beta[i];
  }
  for (int i = 0; i < 3; i++)
    result += beta[i] * u[i];
  return result / sumbeta;
}
inline real weno5_Z(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

  real sumbeta = 0, result = 0;
  real C = 1, qq = 2, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    beta[i] = gamma[i] * (C + pow(tau / (beta[i] + eps), qq));
    sumbeta += beta[i];
  }
  for (int i = 0; i < 3; i++)
    result += beta[i] * u[i];
  return result / sumbeta;
}

inline real Teno5_ZCT4(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau / (beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }
  real CT = 1e-4 * sumbeta;
  // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
  unsigned short flag = 0;
  if (beta[0] < CT)
    flag += 1;
  if (beta[1] < CT)
    flag += 2;
  if (beta[2] < CT)
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}

inline real Teno5_ZCT7(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau / (beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }
  real CT = 1e-7 * sumbeta;
  // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
  unsigned short flag = 0;
  if (beta[0] < CT)
    flag += 1;
  if (beta[1] < CT)
    flag += 2;
  if (beta[2] < CT)
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}

// 计算η值，用于TENO算法
inline double calculate_eta(double delta_f_im12, double delta_f_ip12, double delta_f_ip32) {
    const double xi = 1e-3;       // 小量防止除零
    const double Cr = 0.24;       // 常数参数
    const double epsilon = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;  // ε计算

    // 计算当前点和下一个点的η值
    double eta_i = (std::abs(2 * delta_f_ip12 * delta_f_im12) + epsilon) / 
                  (delta_f_ip12 * delta_f_ip12 + delta_f_im12 * delta_f_im12 + epsilon);
    
    double eta_ip1 = (std::abs(2 * delta_f_ip32 * delta_f_ip12) + epsilon) / 
                    (delta_f_ip32 * delta_f_ip32 + delta_f_ip12 * delta_f_ip12 + epsilon);

    // 返回最小的η值
    return std::min({eta_i, eta_ip1});
}

// TENO-A算法核心函数，返回三个布尔值表示是否使用对应模板
inline std::array<bool, 3> tenoA(double beta0, double beta1, double beta2, double tau, double eta_min) {
    const double epsilon = 1e-40;  // 防止除零的小量
    const double q = 6;            // 指数参数
    const double C = 1;            // 常数
    const double Cr = 0.24;        // 与η计算相同的常数
    const double alpha1 = 10.0;    // 参数1
    const double alpha2 = 5.0;     // 参数2

    // 计算γ值
    double gamma0 = pow(C + tau / (beta0 + epsilon), q);
    double gamma1 = pow(C + tau / (beta1 + epsilon), q);
    double gamma2 = pow(C + tau / (beta2 + epsilon), q);
    double gamma_sum = gamma0 + gamma1 + gamma2;

    // 计算χ值
    double chi0 = gamma0 / gamma_sum;
    double chi1 = gamma1 / gamma_sum;
    double chi2 = gamma2 / gamma_sum;

    // 计算m参数
    double m = 1.0 - std::min(1.0, eta_min / Cr);

    // 计算g(m)和β
    double g_m = pow(1 - m, 4) * (1 + 4 * m);
    double beta = alpha1 - alpha2 * (1 - g_m);

    // 计算阈值CT
    double CT = pow(10, -beta);

    // 返回三个布尔值，表示是否使用对应模板
    return {chi0 > CT, chi1 > CT, chi2 > CT};
}

// 五阶TENO重构函数
inline real Teno5_CongAA(std::array<real, 5> q)
{
     real eps = 1e-40; // 小量防止除零
    
    // 计算三个模板的光滑指示器β
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),
        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),
        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };
    
    real tau = std::abs(beta[2] - beta[0]);  // 计算τ值
    
    // 计算η值 (lambda替代calculate_eta)
    auto calc_eta = [](real a, real b, real c) {
        const real xi = 1e-3, Cr = 0.24;
        const real epsilon = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;
        real eta1 = (std::abs(2*b*a) + epsilon) / (b*b + a*a + epsilon);
        real eta2 = (std::abs(2*c*b) + epsilon) / (c*c + b*b + epsilon);
        return std::min(eta1, eta2);
    };
    
    real eta_min = std::min({
        calc_eta(q[0], q[1], q[2]),
        calc_eta(q[1], q[2], q[3]),
        calc_eta(q[2], q[3], q[4])
    });

    // 内联tenoA逻辑
    const real q_val = 6, C = 1, Cr = 0.24, alpha1 = 10.0, alpha2 = 5.0;
    auto gamma = [&](real b) { return pow(C + tau/(b + eps), q_val); };
    real gamma_sum = gamma(beta[0]) + gamma(beta[1]) + gamma(beta[2]);
    
    real m = 1.0 - std::min(1.0, eta_min / Cr);
    real g_m = pow(1 - m, 4) * (1 + 4 * m);
    real CT = pow(10, -std::round(alpha1 - alpha2*(1 - g_m)));

    // // 直接计算标志位
    // int flag = 0;
    // if (gamma(beta[0])/gamma_sum < CT) flag += 1;
    // if (gamma(beta[1])/gamma_sum < CT) flag += 2;
    // if (gamma(beta[2])/gamma_sum < CT) flag += 4;


    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    // constexpr real CT = std::pow(1.5 * 1e-5, 1.0 / 6.0);
    real CT_1 = 1 - CT;
    // real tau = std::abs(beta[2] - beta[0]);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];


    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;

    // 根据标志位选择重构模板
    switch (flag) {
    case 0:  // 111: 使用所有三个模板
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011: 使用后两个模板
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101: 使用第一和第三个模板
    case 3:  // 001: 只使用第三个模板
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110: 使用前两个模板
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010: 只使用第二个模板
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100: 只使用第一个模板
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 使用中心值
        return q[2];
    }
}
// 【王鸿飞】待研究
inline real Teno5_Z(std::array<real, 5> q) {
  //chenyuqing:最后用的非线性插值函数TENO
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau / (beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }
  real CT = 1e-5 * sumbeta;
  // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);

  unsigned short flag = 0;
  if (beta[0] < CT)
    flag += 1;
  if (beta[1] < CT)
    flag += 2;
  if (beta[2] < CT)
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}
inline real Teno5_ZConvex(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  real sumbeta = 0, result = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau / (beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }

  real CT = 1e-5 * sumbeta, sumGamma = 0;
  for (int i = 0; i < 3; i++) {
    if (beta[i] > CT) {
      // result+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
      result += gamma[i] * u[i];
      sumGamma += gamma[i];
    }
  }
  result /= sumGamma;
  return result;
}

const static real CTi = pow(1.5 * 1e-10, 1.0 / 6.0);
inline real Teno5_CongZ(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};












  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  // inline real CT=0.23050581003334941;//4
  constexpr real CT = 0.15704178024750198; // 5
  // inline real CT=0.08;//5
  // inline real CT=0.10699131939336631;//6
  // inline real CT=0.072892337360747711; //7
  // inline real CT=0.033833625914958219;//9
  // inline real CT=0.023050581003334944;//10
  constexpr real CT_1 = 1 - CT;
  real tau = std::abs(beta[2] -
                      beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
  real rr = CT * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);




  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}


inline real Teno5_CongZCT4(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  constexpr real CT = 0.23050581003334941; // 4
  // inline real CT=0.15704178024750198;//5
  // inline real CT=0.157;//5
  // inline real CT=0.10699131939336631;//6
  // inline real CT=0.072892337360747711; //7
  // inline real CT=0.033833625914958219;//9
  // inline real CT=0.023050581003334944;//10
  constexpr real CT_1 = 1 - CT;
  real tau = std::abs(beta[2] -
                      beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
  real rr = CT * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}
inline real Teno5_CongZCT7(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  // inline real CT=0.23050581003334941;//4
  // inline real CT=0.15704178024750198;//5
  // inline real CT=0.157;//5
  // inline real CT=0.10699131939336631;//6
  constexpr real CT = 0.072892337360747711; // 7
  // inline real CT=0.033833625914958219;//9
  //  inline real CT=0.023050581003334944;//10
  constexpr real CT_1 = 1 - CT;
  real tau = std::abs(beta[2] -
                      beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
  real rr = CT * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}

inline real Teno5_CongA(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  // std::array<real,3> beta={
  //         pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2),
  //         pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2),
  //         pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  // real tau=std::abs(beta[2]+beta[0]-2*beta[1]);
  real tau = std::abs(beta[2] - beta[0]);

  // inline real CT=0.23050581003334941;//4
  // inline real CT=0.15704178024750198;//5
  constexpr real CT = 0.1; // 6
  // inline real CT=0.072892337360747711;//7
  // inline real CT=0.033833625914958219;//9
  //  inline real CT=0.023050581003334944;//10
  //  real CT=0.15704178024750198;//5
  real CT_1 = 1 - CT;

  // real
  // rr=CT/1/(beta[minBeta]+eps)-CT_1/(tau+eps);//CT*tau-CT_1*beta[minBeta];
  real mulbeta = beta[0] * beta[1] * beta[2];
  real rr =
      CT * tau * (beta[0] * beta[1] + beta[1] * beta[2] + beta[0] * beta[2]) -
      CT_1 * mulbeta;      // CT*tau-CT_1*beta[minBeta];
  real ll = tau * mulbeta; // tau*beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}
inline real Teno5_CongC(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  // real tau=std::abs(beta[2]+beta[0]-2*beta[1]);

  real tau = std::abs(beta[2] - beta[0]);

  // inline real CT=0.23050581003334941;//4
  // inline real CT=0.15704178024750198;//5
  // inline real CT=0.10699131939336631;//6
  // inline real CT=0.072892337360747711; //7
  // inline real CT=0.033833625914958219;//9
  //  inline real CT=0.023050581003334944;//10
  real CT = 0.15704178024750198; // 5
  real CT_1 = 1 - CT;

  real rr = CT * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];

  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

  // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

  real sumbeta = 0, result = 0, sumGamma = 0;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  for (int i = 0; i < 3; i++) {
    if (ll >= rr * beta[i]) {
      sumGamma += gamma[i];
      result += gamma[i] * u[i];
    }
  }
  return result / sumGamma;
}

inline real Teno5_Cong(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0, result = 0;

  real minBeta = *std::min_element(beta.begin(), beta.end());

  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};
  real CT = 3.5, sumGamma = 0;
  for (int i = 0; i < 3; i++) {
    if (beta[i] < CT * (minBeta + eps)) {
      sumGamma += gamma[i];
      result += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
    }
  }
  result /= sumGamma; // 除法赋值 即 result = result/sumGamma;

  return result;
}

inline real Teno5_Cong2(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0, result = 0;

  real minBeta = *std::min_element(beta.begin(), beta.end());
  if (minBeta < eps)
    return q[2];
  real CT = 2, sumGamma = 0;
  int flag = 0, flags[3] = {1, 2, 4};

  for (int i = 0; i < 3; i++) {
    if (beta[i] < CT * (minBeta + eps)) {
      flag += flags[i];
    }
  }
  switch (flag) {
  case 1:
    /* 1,0,0 */
    result = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  case 2:
    /* 0,1,0 */
    result = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 3:
    /* 1,1,0 */
    result = 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
             5.0 / 16.0 * q[3];
    break;
  case 4:
  case 5:
    /* 0,0,1 */
    /* 1,0,1 */
    result = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 6:
    /* 0,1,1 */
    result = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
             1.0 / 16.0 * q[4];
    break;
  case 7:
    /* 1,1,1 */
    result = 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
             15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  default:
    /* 0,0,0 */
    result = q[2];
    break;
  }

  return result;
}

inline std::array<real, 3> Teno5_BVDCong(std::array<real, 5> q, bool &flag) {
  real eps = 1e-40;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  std::array<real, 3> result = {0, 0};
  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    beta[i] = (C + pow(tau / (beta[i] + eps), qq));
    sumbeta += beta[i];
  }
  for (int i = 0; i < 3; i++)
    beta[i] /= sumbeta;

  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};
  real CT1 = 1e-7, CT2 = 1e-3, sumGamma = 0;
  for (int i = 0; i < 3; i++) {
    if (beta[i] > CT1) {
      sumGamma += gamma[i];
      result[0] += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
    }
  }
  result[0] /= sumGamma;

  real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
  if (!(extremPoint >= 1 || extremPoint <= 0)) {
    result[1] = q[2];
    result[2] = result[0];
    flag = false;
  } else {
    flag = true;
    sumGamma = 0;
    for (int i = 0; i < 3; i++) {
      if (beta[i] > CT2) {
        sumGamma += gamma[i];
        result[1] += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
      }
    }
    result[1] /= sumGamma;
    // auto minBetai=std::max_element(beta.begin(),beta.end())-beta.begin();
    // result[1]=(*u[minBetai])(q[minBetai],q[minBetai+1],q[minBetai+2]);
    // {result[2]=q[2];}
  }

  return result;
}

inline std::array<real, 3> Teno5_BVDMR(std::array<real, 5> q, bool &flag) {
  real eps = 1e-6;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  std::array<real, 3> result = {0, 0};
  real sumbeta = 0;

  auto minBetap = std::min_element(beta.begin(), beta.end());
  auto minBetai = minBetap - beta.begin();
  auto minBeta = *minBetap;
  if (minBeta < eps)
    return {q[2], q[2], q[2]};

  real CT1 = 20, CT2 = 1.5, sumGamma = 0;
  int flag1 = 0, flag2 = 0, flags[] = {1, 2, 4};
  for (int i = 0; i < 3; i++) {
    if (beta[i] < CT1 * (minBeta + eps)) {
      flag1 += flags[i];
    }
    if (beta[i] <
        CT2 * (minBeta +
               eps)) //(extremPoint>=2-minBetai || extremPoint<=1-minBetai))
    {
      flag2 += flags[i];
    }
  }
  switch (flag1) {
  case 1:
    /* 1,0,0 */
    result[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  case 2:
    /* 0,1,0 */
    result[0] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 3:
    /* 1,1,0 */
    result[0] = 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
                5.0 / 16.0 * q[3];
    break;
  case 4:
  case 5:
    /* 0,0,1 */
    /* 1,0,1 */
    result[0] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 6:
    /* 0,1,1 */
    result[0] = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
                1.0 / 16.0 * q[4];
    break;
  case 7:
    /* 1,1,1 */
    result[0] = 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
                15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  default:
    /* 0,0,0 */
    result[0] = q[2];
    break;
  }
  switch (flag2) {
  case 1:
    /* 1,0,0 */
    result[1] = -5.0 / 12.0 * q[0] + 1.0 / 3.0 * q[1] + 13.0 / 12.0 * q[2];
    break;
  case 2:
    /* 0,1,0 */
    result[1] = 1.0 / 12.0 * q[1] + 1.0 / 3.0 * q[2] + 7.0 / 12.0 * q[3];
    break;
  case 3:
    /* 1,1,0 */
    result[1] = -9.0 / 80.0 * q[0] - 17.0 / 80.0 * q[1] + 33.0 / 80.0 * q[2] +
                39.0 / 80.0 * q[3];
    break;
  case 4:
  case 5:
    /* 0,0,1 */
    /* 1,0,1 */
    result[1] = 7.0 / 12.0 * q[2] + 1.0 / 3.0 * q[3] + 1.0 / 12.0 * q[4];
    break;
  case 6:
    /* 0,1,1 */
    result[1] = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
                1.0 / 16.0 * q[4];
    break;
  case 7:
    /* 1,1,1 */
    result[1] = -3.0 / 160.0 * q[0] + 1.0 / 80.0 * q[1] + 9.0 / 20.0 * q[2] +
                51.0 / 80.0 * q[3] - 13.0 / 160.0 * q[4];
    break;
  default:
    /* 0,0,0 */
    result[1] = q[2];
    break;
  }
  result[1] = result[0];

  return result;
}

inline real Teno5_CongSort(std::array<real, 5> q) {
  real eps = 1e-12;
  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};

  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real CC = (beta[1] * 2 - beta[2] - beta[0]) / 2;
  // beta[0]+=2.0/3.0*CC;
  // beta[1]-=1.0/3.0*CC;
  // beta[2]+=2.0/3.0*CC;

  // 排序
  std::array<real, 3> index = {0, 1, 2};
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (beta[a] < beta[b]); });
  // if(beta[index[2]]<eps) return q[2];
  int ii = index[0];
  real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];

  real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
  real extremPoint2 = 0.5 + (q[2] - q[3]) / (q[2] - 2 * q[3] + q[4] + eps);

  bool flags[3] = {
      true //(u1<*std::max_element(q.begin(),q.end())&&(u1>*std::min_element(q.begin(),q.end())))
      ,
      (extremPoint >= 1 || extremPoint <= 0),
      (extremPoint2 >= 1 || extremPoint2 <= 0)}; //;|| (factor>6);
  // bool flag=flags[1]||flags[2];//;|| (factor>6);
  real CT2, CT = 3, sumGamma = gamma[index[0]];
  // if(!flag)
  // {
  //     CT2=3;
  //     CT=CT2;//-beta[index[1]]/(beta[index[0]]+eps)/1.5;
  // }
  // else
  // {
  //     CT2=3;
  //     CT=CT2;//-beta[index[1]]/(beta[index[0]]+eps)/1.5;
  // }
  // return result/sumGamma;

  ii = index[1];
  if (beta[ii] < CT * (beta[index[0]] + eps)) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  } else {
    // real extremPoint=0.5*(q[1]-q[3])/(q[1]-2*q[2]+q[3]+eps);
    //  int i0=index[0];
    //  if(!flags[i0]) return q[2];
    return result / sumGamma;
  }
  // {if(flag) return result/sumGamma;
  //  else
  //  return q[2];}

  ii = index[2];
  if (beta[ii] < (CT) * (beta[index[0]] + eps)) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  }
  result /= sumGamma;

  return result;
}

inline real Teno5_CongSortPositive(std::array<real, 5> q) {
  real eps = 1e-12;
  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};

  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  // 排序
  std::array<real, 3> index = {0, 1, 2};
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (beta[a] < beta[b]); });
  // if(beta[index[2]]<eps) return q[2];
  int ii = index[0];
  real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];
  if (result < 0)
    return q[2];

  real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
  real extremPoint2 = 0.5 + (q[2] - q[3]) / (q[2] - 2 * q[3] + q[4] + eps);
  real factor = std::abs((q[2] - q[1] + eps) / (q[2] - q[3] + eps));
  bool flag =
      (extremPoint >= 1 + eps && extremPoint >= 0 - eps) ||
      (extremPoint2 <= 1 + eps && extremPoint2 >= 0 - eps); //;|| (factor>6);
  real CT2 = factor, CT, sumGamma = gamma[index[0]];
  if (flag) {
    CT2 = 4;
    CT = CT2 - beta[index[1]] / (beta[index[0]] + eps);
  } else {
    CT2 = 20;
    CT = CT2 - beta[index[1]] / (beta[index[0]] + eps) / 1.5;
  }

  // return result/sumGamma;

  ii = index[1];
  real u1 = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  if (beta[ii] < CT * (beta[index[0]]) && u1 > 0) {
    sumGamma += gamma[ii];
    result += gamma[ii] * u1;
  } else {
    if (!flag)
      return result / sumGamma;
    else
      return q[2];
  }

  ii = index[2];
  real u2 = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  if (beta[ii] < (CT) * (beta[index[0]] + eps) && u1 > 0) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  }
  result /= sumGamma;

  return result;
}

inline real Teno5_CongSortabs(std::array<real, 5> q) {
  real eps = 1e-12;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  // 排序
  std::array<real, 3> index = {0, 1, 2};
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (beta[a] < beta[b]); });

  std::array<real (*)(real, real, real), 3> u = {&u1, &u2, &u3};

  real factor = std::abs((q[2] - q[1] + eps) / (q[2] - q[3] + eps));
  real CT2 = factor, CT, sumGamma = gamma[index[0]];
  if (factor > 6) {
    CT2 = 6;
    CT = CT2 - beta[index[1]] / (beta[index[0]] + eps);
  } else {
    CT2 = 20;
    CT = CT2 - beta[index[1]] / (beta[index[0]] + eps) / 1.5;
  }
  int ii = index[0];
  real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];
  // return result/sumGamma;

  if (beta[index[2]] < eps)
    return q[2];
  // real critical=pow(q[2]-q[1],2);
  // CT2=std::max(20.0*pow(critical,2.0),4.0);

  // CT=4;
  ii = index[1];
  if (beta[ii] < CT * (beta[index[0]])) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  } else {
    if (factor < 10)
      return result / sumGamma;
    else
      return q[2];
  }

  ii = index[2];
  if (beta[ii] < (CT) * (beta[index[0]] + eps)) {
    sumGamma += gamma[ii];
    result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
  }
  result /= sumGamma;

  return result;
}

inline real Teno5_CongIncrease(std::array<real, 5> q) {
  real eps = 1e-20;
  std::array<real, 3> gamma = {1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0};
  std::array<real, 3> beta;
  beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  // 排序
  std::array<real, 3> index = {0, 1, 2};
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (beta[a] < beta[b]); });

  // 过于光滑
  if (beta[index[2]] <= 1e-10)
    return q[2];

  std::array<real, 3> u;
  u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
  u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
  real u0 = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
  real sumbeta = beta[index[0]], sumGamma = gamma[index[0]];
  real result = u[index[0]] * gamma[index[0]];

  int ii = index[1];
  real CT = 10, a1 = (beta[index[0]] + eps);
  if (beta[ii] < CT * a1) {
    sumGamma += gamma[ii];
    result += gamma[ii] * u[ii];
    sumbeta += beta[ii];
  } else {
    return result / sumGamma;
  }

  real CT2 = beta[ii];
  ii = index[2];
  real ddd = 3;
  if (pow(beta[ii] * a1, ddd) - 5 * pow(CT2 * CT2, ddd) < 0) {
    sumGamma += gamma[ii];
    result += gamma[ii] * u[ii];
    sumbeta += beta[ii];
  } else {
    return result / sumGamma;
  }

  return result / sumGamma;
}


inline real musclInterpolation(real q1, real q2, real q3) {

  real delta;
  real deltam, deltap;
  deltam = q2 - q1;
  deltap = q3 - q2;

  // minmod
  real beta = 1.0;
  if (deltap > 0) {
    delta = std::max(0.0, std::max(std::min(beta * deltam, deltap),
                                   std::min(deltam, beta * deltap)));
  } else {
    delta = std::min(0.0, std::min(std::max(beta * deltam, deltap),
                                   std::max(deltam, beta * deltap)));
  }
  return q2 + delta * 0.5;
}
inline real musclIn5(std::array<real, 5> q) {
  return musclInterpolation(q[1],q[2],q[3]);
}

inline std::array<real, 2> THINC(real q1, real q2, real q3) {
  if ((q1 - q2) * (q2 - q3) < 0)
    return {q2, q2};

  real qmax, qmin;
  if (q1 > q3) {
    qmax = q1;
    qmin = q3;
  } else {
    qmax = q3;
    qmin = q1;
  }
  real beta = 2;
  real T1 = exp(2.0 * beta), T3 = exp(-2.0 * beta);
  real f = (qmax - q2) / (q2 - qmin);
  T1 *= f;
  T3 *= f;
  real qBar = (q1 + q3) / 2, dq = (q1 - q3) / 2;
  return {qBar + dq * ((1 - T1) / (1 + T1)), qBar + dq * ((1 - T3) / (1 + T3))};
}

inline real THINC1(real q1, real q2, real q3) {
  if ((q1 - q2) * (q2 - q3) < 1e-20)
    return q2;

  real qmax, qmin;
  if (q1 > q3) {
    qmax = q1;
    qmin = q3;
  } else {
    qmax = q3;
    qmin = q1;
  }
  real beta = 1.2;
  real T1 = exp(2.0 * beta), T3 = exp(-2.0 * beta);
  real f = (qmax - q2) / (q2 - qmin);
  T1 *= f;
  T3 *= f;
  real qBar = (q1 + q3) / 2, dq = (q1 - q3) / 2;
  return qBar + dq * ((1 - T3) / (1 + T3));
}

inline real Linear5(std::array<real, 5> q) {
  return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
         15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
}



//【王鸿飞】begin插值格式开发（1.0）手搓
inline real whf_TcnsN(std::array<real, 5> q)
{
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
              1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 
              1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 
              1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);

    // 求光滑度量（带入光滑因子）
    const real C = 1.0;
    const real eps = 1e-40;
    const int q_val = 6;
    auto gamma_cal = [&](real a) { return pow(C + tau/(a + eps), q_val); };

    std::array<real, 3> gamma;

    gamma[0] = gamma_cal(beta[0]);
    gamma[1] = gamma_cal(beta[1]);
    gamma[2] = gamma_cal(beta[2]);

    // 光滑度量，归一化
    real gamma_sum = gamma[0] + gamma[1] + gamma[2];

    std::array<real, 3> chi;
    chi[0] = gamma[0]/gamma_sum;
    chi[1] = gamma[1]/gamma_sum;
    chi[2] = gamma[2]/gamma_sum;

    // 用光滑度量，构造截止函数
    real CT = 1e-10;

    unsigned short flag = 0;
    if (chi[0] < CT) flag += 1;
    if (chi[1] < CT) flag += 2;
    if (chi[2] < CT) flag += 4;

    // 插值
    switch (flag) {
    case 0:  // 111
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101
    case 3:  // 001
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}

inline real whf_TcnsN_S(std::array<real, 5> q)
{
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
              1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 
              1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 
              1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);

    // 辅助间断函数
    real CT = 1e-5;

    const real q_A = 6;
    real CT_p = pow(1.5*CT/(1.0 - CT) , 1.0/q_A);

    // 用辅助间断函数，构造截止函数
    real minBeta = std::min({beta[0], beta[1], beta[2]});
    // real minBeta = std::min({beta[0], beta[1], beta[2]});
    real CC = 1.0;
    
    std::array<real, 3> leftleft;
    std::array<real, 3> rightright;

    leftleft[0] = beta[0]*tau;
    leftleft[1] = beta[1]*tau;
    leftleft[2] = beta[2]*tau;

    rightright[0] = (CT_p*tau - CC*(1.0 - CT_p)*beta[0])*minBeta;
    rightright[1] = (CT_p*tau - CC*(1.0 - CT_p)*beta[1])*minBeta;
    rightright[2] = (CT_p*tau - CC*(1.0 - CT_p)*beta[2])*minBeta;

    unsigned short flag = 0;
    if (leftleft[0] < rightright[0]) flag += 1;
    if (leftleft[1] < rightright[1]) flag += 2;
    if (leftleft[2] < rightright[2]) flag += 4;

    // 插值
    switch (flag) {
    case 0:  // 111
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101
    case 3:  // 001
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}



inline real whf_TcnsN_AS(std::array<real, 5> q)
{
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0/1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
              1.0/4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0/1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 
              1.0/4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0/1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 
              1.0/4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);

    // 求光滑度量（带入光滑因子）
    const real C = 1.0;
    const real eps = 1e-40;
    const real q_val = 6;
    auto gamma_cal = [&](real a) { return pow(C + tau/(a + eps), q_val); };

    std::array<real, 3> gamma;

    gamma[0] = gamma_cal(beta[0]);
    gamma[1] = gamma_cal(beta[1]);
    gamma[2] = gamma_cal(beta[2]);

    // 光滑度量，归一化
    real gamma_sum = gamma[0] + gamma[1] + gamma[2];

    std::array<real, 3> chi;
    chi[0] = gamma[0]/gamma_sum;
    chi[1] = gamma[1]/gamma_sum;
    chi[2] = gamma[2]/gamma_sum;

    // 计算自适应阈值CT_A
    const real C_r = 0.24, alpha1 = 10.0, alpha2 = 5.0;  //参数

    auto eta_calc = [](real a, real b, real c) { // 引入比率函数 η_k
        //calc为calculate的缩写
        const real xi = 1e-3, C_r_inner = 0.24;
        const real epsilon_A = (0.9*C_r_inner)/(1.0 - 0.9*C_r_inner)*xi*xi;
        real eta1 = (std::abs(2*b*a) + epsilon_A)/(b*b + a*a + epsilon_A);
        real eta2 = (std::abs(2*c*b) + epsilon_A)/(c*c + b*b + epsilon_A);
        return std::min(eta1, eta2);
    };
    
    real eta_min = std::min({
        eta_calc(q[0], q[1], q[2]),
        eta_calc(q[1], q[2], q[3]),
        eta_calc(q[2], q[3], q[4])
    });     // 计算η_min值
    real m = 1.0 - std::min(1.0, eta_min/C_r);

    real g_m = pow(1 - m, 4)*(1 + 4*m);
    int beta_A = std::floor(alpha1 - alpha2*(1 - g_m));
    real CT_A = pow(10, -beta_A);

    // 辅助间断函数
    const int q_A = 6;
    real CT_p = pow(1.5*CT_A/(1.0 - CT_A) , 1.0/q_A);

    // 用辅助间断函数，构造截止函数
    real minBeta = std::min({beta[0], beta[1], beta[2]});
    real CC = 1.0;
    
    std::array<real, 3> leftleft;
    std::array<real, 3> rightright;

    leftleft[0] = beta[0]*tau;
    leftleft[1] = beta[1]*tau;
    leftleft[2] = beta[2]*tau;

    rightright[0] = (CT_p*tau - CC*(1.0 - CT_p)*beta[0])*minBeta;
    rightright[1] = (CT_p*tau - CC*(1.0 - CT_p)*beta[1])*minBeta;
    rightright[2] = (CT_p*tau - CC*(1.0 - CT_p)*beta[2])*minBeta;

    unsigned short flag = 0;
    if (leftleft[0] < rightright[0]) flag += 1;
    if (leftleft[1] < rightright[1]) flag += 2;
    if (leftleft[2] < rightright[2]) flag += 4;

    // 插值
    switch (flag) {
    case 0:  // 111
        return 3.0/128.0 * q[0] - 5.0/32.0 * q[1] + 45.0/64.0 * q[2] + 15.0/32.0 * q[3] - 5.0/128.0 * q[4];
    case 1:  // 011
        return -1.0/16.0 * q[1] + 9.0/16.0 * q[2] + 9.0/16.0 * q[3] - 1.0/16.0 * q[4];
    case 2:  // 101
    case 3:  // 001
        return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    case 4:  // 110
        return 1.0/16.0 * q[0] - 5.0/16.0 * q[1] + 15.0/16.0 * q[2] + 5.0/16.0 * q[3];
    case 5:  // 010
        return -1.0/8.0 * q[1] + 3.0/4.0 * q[2] + 3.0/8.0 * q[3];
    case 6:  // 100
        return 3.0/8.0 * q[0] - 5.0/4.0 * q[1] + 15.0/8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}

inline real whf_TcnsN_LAD(std::array<real, 5> q)
{
    // （步骤1）光滑因子β
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);
    // （步骤2）计算LAD自适应阈值CT
    real eps_lad = 1e-40;
    std::array<real, 3> chi;
    chi[0] = tau / (beta[0] + eps_lad);
    chi[1] = tau / (beta[1] + eps_lad);
    chi[2] = tau / (beta[2] + eps_lad);

    real chi_max = std::max({chi[0], chi[1], chi[2]});

    real hh = 10.0;
    real theta = 1.0 / (1.0 + chi_max/hh);
    int b_l = 4;
    int b_u = 10;
    // 预计算CT值表
    const std::array<real, 7> CT_table = {  // m从4到10对应的CT值
        1e-4,   // m=4
        1e-5,   // m=5 
        1e-6,   // m=6
        1e-7,   // m=7
        1e-8,   // m=8
        1e-9,   // m=9
        1e-10   // m=10
    };
    int m = b_l + std::floor(theta*(b_u - b_l));
    m = std::clamp(m, b_l, b_u);  // 确保m在4-10范围内
    real CT = CT_table[m - b_l];   // 直接查表获取CT值

    // （步骤3）截断函数
    // 内联tenoA逻辑
    const real C = 1.0;
    const real eps = 1e-40;
    const int q_val = 6;

    auto gamma_cal = [&](real a) { return pow(C + tau/(a + eps), q_val); };

    std::array<real, 3> gamma;
    gamma[0] = gamma_cal(beta[0]);
    gamma[1] = gamma_cal(beta[1]);
    gamma[2] = gamma_cal(beta[2]);

    real gamma_sum = gamma[0] + gamma[1] + gamma[2];
    // （步骤4）计算模板选择标志位
    unsigned short flag = 0;
    if (gamma[0]/gamma_sum < CT) flag += 1;
    if (gamma[1]/gamma_sum < CT) flag += 2;
    if (gamma[2]/gamma_sum < CT) flag += 4;

    // （步骤5）基于标志位选择重构模板（执行插值计算）
    switch (flag) {
    case 0:  // 111: 全模板(5阶)
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:  // 011: 右侧2模板(4阶)
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:  // 101: 不连续区域 -> 使用case3
    case 3:  // 001: 仅模板2(3阶)
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:  // 110: 左侧2模板(4阶)
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:  // 010: 仅模板1(3阶)
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    case 6:  // 100: 仅模板0(3阶)
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}


//【王鸿飞】end插值格式开发（1.0）手搓


//【王鸿飞】begin插值格式开发（2.0）AI-原格式
inline real congTcns5_ZCT5(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0/1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0/4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0/1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0/4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0/1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0/4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau/(beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }
  real CT = 1e-5 * sumbeta;
  // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
  unsigned short flag = 0;
  if (beta[0] < CT)
    flag += 1;
  if (beta[1] < CT)
    flag += 2;
  if (beta[2] < CT)
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0/128.0 * q[0] - 5.0/32.0 * q[1] + 45.0/64.0 * q[2] +
           15.0/32.0 * q[3] - 5.0/128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0/16.0 * q[1] + 9.0/16.0 * q[2] + 9.0/16.0 * q[3] -
           1.0/16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0/16.0 * q[0] - 5.0/16.0 * q[1] + 15.0/16.0 * q[2] +
           5.0/16.0 * q[3];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0/8.0 * q[1] + 3.0/4.0 * q[2] + 3.0/8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0/8.0 * q[0] - 5.0/4.0 * q[1] + 15.0/8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}

inline real congTcns5_ZCT10(std::array<real, 5> q) {
  real eps = 1e-40;
  std::array<real, 3> beta;
  beta[0] = 1.0/1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0/4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

  beta[1] = 1.0/1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0/4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

  beta[2] = 1.0/1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0/4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

  real sumbeta = 0;
  real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
  for (int i = 0; i < 3; i++) {
    real tempp = C + tau/(beta[i] + eps);
    tempp *= tempp;
    beta[i] = tempp * tempp * tempp;
    sumbeta += beta[i];
  }
  real CT = 1e-10 * sumbeta;
  // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
  unsigned short flag = 0;
  if (beta[0] < CT)
    flag += 1;
  if (beta[1] < CT)
    flag += 2;
  if (beta[2] < CT)
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0/128.0 * q[0] - 5.0/32.0 * q[1] + 45.0/64.0 * q[2] +
           15.0/32.0 * q[3] - 5.0/128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0/16.0 * q[1] + 9.0/16.0 * q[2] + 9.0/16.0 * q[3] -
           1.0/16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0/16.0 * q[0] - 5.0/16.0 * q[1] + 15.0/16.0 * q[2] +
           5.0/16.0 * q[3];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0/8.0 * q[1] + 3.0/4.0 * q[2] + 3.0/8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0/8.0 * q[0] - 5.0/4.0 * q[1] + 15.0/8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}



inline real whf_ai_TcnsN_S_1(std::array<real, 5> q) {
    // 常量定义
    const real C_T = 1e-5;
    const real q_A = 6.0;
    const real C = 1.0;
    
    // 计算辅助截断函数 C_T_prime
    const real C_T_prime = std::pow(1.5 * C_T / (1.0 - C_T), 1.0 / q_A);
    
    // 解包输入值
    const real q_im2 = q[0]; // q_{i-2}
    const real q_im1 = q[1]; // q_{i-1}
    const real q_i   = q[2]; // q_i
    const real q_ip1 = q[3]; // q_{i+1}
    const real q_ip2 = q[4]; // q_{i+2}
    
    // 计算局部光滑因子 β
    const real beta0 = 1.0 / 1.0 * std::pow(1.0 * q_im2 - 2.0 * q_im1 + 1.0 * q_i, 2) +
              1.0 / 4.0 * std::pow(1.0 * q_im2 - 4.0 * q_im1 + 3.0 * q_i, 2);
    const real beta1 = 1.0 / 1.0 * std::pow(1.0 * q_im1 - 2.0 * q_i + 1.0 * q_ip1, 2) + 
              1.0 / 4.0 * std::pow(1.0 * q_im1 + 0.0 * q_i - 1.0 * q_ip1, 2);
    const real beta2 = 1.0 / 1.0 * std::pow(1.0 * q_i - 2.0 * q_ip1 + 1.0 * q_ip2, 2) + 
              1.0 / 4.0 * std::pow(3.0 * q_i - 4.0 * q_ip1 + 1.0 * q_ip2, 2);
        
    // 计算全局光滑因子 τ
    const real tau = std::fabs(beta2 - beta0);
    
    // 计算最小光滑因子 β_j
    // const real beta_j = std::min(std::min(beta0, beta1), beta2);
    const real beta_j = std::min({beta0, beta1, beta2});
    
    // 计算δ_k
    int delta0, delta1, delta2;
    
    // δ0计算
    if (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta0) {
        delta0 = 0;
    } else {
        delta0 = 1;
    }
    
    // δ1计算
    if (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta1) {
        delta1 = 0;
    } else {
        delta1 = 1;
    }
    
    // δ2计算
    if (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta2) {
        delta2 = 0;
    } else {
        delta2 = 1;
    }
    
    // 系数表 [a_{i-2}, a_{i-1}, a_i, a_{i+1}, a_{i+2}]
    const std::array<std::array<real, 5>, 8> coeff_table = {{
        // δ0, δ1, δ2 = (1,1,1)
        { 3.0/128, -5.0/32, 45.0/64, 15.0/32, -5.0/128 },
        // (0,1,1)
        { 0.0, -1.0/16, 9.0/16, 9.0/16, -1.0/16 },
        // (1,0,1)
        { 0.0, 0.0, 3.0/8, 3.0/4, -1.0/8 },
        // (0,0,1)
        { 0.0, 0.0, 3.0/8, 3.0/4, -1.0/8 },
        // (1,1,0)
        { 1.0/16, -5.0/16, 15.0/16, 5.0/16, 0.0 },
        // (0,1,0)
        { 0.0, -1.0/8, 3.0/4, 3.0/8, 0.0 },
        // (1,0,0)
        { 3.0/8, -5.0/4, 15.0/8, 0.0, 0.0 },
        // (0,0,0)
        { 0.0, 0.0, 1.0, 0.0, 0.0 }
    }};
    
    // 确定系数索引
    size_t index = 0;
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) index = 0;
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) index = 1;
    else if (delta0 == 1 && delta1 == 0 && delta2 == 1) index = 2;
    else if (delta0 == 0 && delta1 == 0 && delta2 == 1) index = 3;
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) index = 4;
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) index = 5;
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) index = 6;
    else index = 7; // (0,0,0)
    
    // 计算重构值
    real f_hat = 0.0;
    for (int i = 0; i < 5; ++i) {
        f_hat += coeff_table[index][i] * q[i];
    }
    
    return f_hat;
}

inline real whf_ai_TcnsN_S_2(std::array<real, 5> q) {
    // 解包输入值
    real q_im2 = q[0]; // q_{i-2}
    real q_im1 = q[1]; // q_{i-1}
    real q_i   = q[2]; // q_i
    real q_ip1 = q[3]; // q_{i+1}
    real q_ip2 = q[4]; // q_{i+2}

    // 计算局部光滑因子
    real beta0 = std::pow(q_im2 - 2*q_im1 + q_i, 2) 
               + 0.25 * std::pow(q_im2 - 4*q_im1 + 3*q_i, 2);
    
    real beta1 = std::pow(q_im1 - 2*q_i + q_ip1, 2) 
               + 0.25 * std::pow(q_im1 - q_ip1, 2);
    
    real beta2 = std::pow(q_i - 2*q_ip1 + q_ip2, 2) 
               + 0.25 * std::pow(3*q_i - 4*q_ip1 + q_ip2, 2);

    // 计算全局光滑因子
    real tau = std::fabs(beta2 - beta0);

    // 计算辅助截断函数
    const real C_T = 1e-5;
    const real q_A = 6.0;
    real C_T_prime = std::pow(1.5 * C_T / (1.0 - C_T), 1.0/q_A);

    // 找到最小光滑因子
    real beta_j = std::min({beta0, beta1, beta2});

    // 计算截断函数
    int delta0 = 0, delta1 = 0, delta2 = 0;
    const real C = 1.0;
    
    // 计算delta_k的条件
    auto condition = [&](real beta_k) {
        return beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta_k;
    };

    delta0 = condition(beta0) ? 0 : 1;
    delta1 = condition(beta1) ? 0 : 1;
    delta2 = condition(beta2) ? 0 : 1;

    // 定义系数表 [a_{i-2}, a_{i-1}, a_i, a_{i+1}, a_{i+2}]
    const std::array<std::array<real, 5>, 8> coefficients = {{
        // δ0, δ1, δ2 = 1,1,1
        { 3.0/128, -5.0/32, 45.0/64, 15.0/32, -5.0/128 },
        // 0,1,1
        { 0, -1.0/16, 9.0/16, 9.0/16, -1.0/16 },
        // 1,0,1
        { 0, 0, 3.0/8, 3.0/4, -1.0/8 },
        // 0,0,1
        { 0, 0, 3.0/8, 3.0/4, -1.0/8 },
        // 1,1,0
        { 1.0/16, -5.0/16, 15.0/16, 5.0/16, 0 },
        // 0,1,0
        { 0, -1.0/8, 3.0/4, 3.0/8, 0 },
        // 1,0,0
        { 3.0/8, -5.0/4, 15.0/8, 0, 0 },
        // 0,0,0
        { 0, 0, 1.0, 0, 0 }
    }};

    // 根据δ值选择系数
    int index = 0;
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) index = 0;
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) index = 1;
    else if (delta0 == 1 && delta1 == 0 && delta2 == 1) index = 2;
    else if (delta0 == 0 && delta1 == 0 && delta2 == 1) index = 3;
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) index = 4;
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) index = 5;
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) index = 6;
    else index = 7;  // 0,0,0

    // 计算数值通量
    const auto& a = coefficients[index];
    real f_hat = a[0] * q_im2 
               + a[1] * q_im1 
               + a[2] * q_i 
               + a[3] * q_ip1 
               + a[4] * q_ip2;

    return f_hat;
}

inline real whf_ai_TcnsN_A_1(std::array<real, 5> q) {
    // 常量定义
    const real C = 1.0;
    const real f_exp = 6.0;  // 指数q
    const real epsilon = 1e-40;
    const real C_r = 0.24;
    const real xi = 1e-3;
    const real epsilonA = (0.9 * C_r) / (1.0 - 0.9 * C_r) * xi * xi;
    
    const real alpha1 = 10.0;  // 
    const real alpha2 = 5.0;  // 

    // 解包输入值
    const real f_im2 = q[0]; // f_{i-2}
    const real f_im1 = q[1]; // f_{i-1}
    const real f_i   = q[2]; // f_i
    const real f_ip1 = q[3]; // f_{i+1}
    const real f_ip2 = q[4]; // f_{i+2}

    // 计算差分
    const real df_im3_2 = f_im1 - f_im2; // Δf_{i-3/2}
    const real df_im1_2 = f_i - f_im1;   // Δf_{i-1/2}
    const real df_ip1_2 = f_ip1 - f_i;   // Δf_{i+1/2}
    const real df_ip3_2 = f_ip2 - f_ip1; // Δf_{i+3/2}

    // 计算光滑度量β
    const real beta0 = std::pow(f_im2 - 2*f_im1 + f_i, 2) + 
                      0.25 * std::pow(f_im2 - 4*f_im1 + 3*f_i, 2);
    
    const real beta1 = std::pow(f_im1 - 2*f_i + f_ip1, 2) + 
                      0.25 * std::pow(f_im1 - f_ip1, 2);
    
    const real beta2 = std::pow(f_i - 2*f_ip1 + f_ip2, 2) + 
                      0.25 * std::pow(3*f_i - 4*f_ip1 + f_ip2, 2);

    // 计算全局光滑因子τ
    const real tau = std::abs(beta2 - beta0);

    // 计算γ_k
    const real gamma0 = std::pow(C + tau/(beta0 + epsilon), f_exp);
    const real gamma1 = std::pow(C + tau/(beta1 + epsilon), f_exp);
    const real gamma2 = std::pow(C + tau/(beta2 + epsilon), f_exp);

    // 归一化光滑度量χ_k
    const real sum_gamma = gamma0 + gamma1 + gamma2;
    const real chi0 = gamma0 / sum_gamma;
    const real chi1 = gamma1 / sum_gamma;
    const real chi2 = gamma2 / sum_gamma;

    // 计算η值
    const real eta_im1 = (std::abs(2*df_im1_2*df_im3_2) + epsilonA) /
                         (std::pow(df_im1_2, 2) + std::pow(df_im3_2, 2) + epsilonA);
    
    const real eta_i = (std::abs(2*df_ip1_2*df_im1_2) + epsilonA) /
                      (std::pow(df_ip1_2, 2) + std::pow(df_im1_2, 2) + epsilonA);
    
    const real eta_ip1 = (std::abs(2*df_ip3_2*df_ip1_2) + epsilonA) /
                         (std::pow(df_ip3_2, 2) + std::pow(df_ip1_2, 2) + epsilonA);

    // 计算m
    const real eta_min = std::min({eta_im1, eta_i, eta_ip1});
    const real m = 1.0 - std::min(1.0, eta_min/C_r);

    // 计算g(m)
    const real g_m = std::pow(1.0 - m, 4) * (1.0 + 4.0*m);

    // 计算β_A和C_TA
    const real beta_A = alpha1 - alpha2 * (1.0 - g_m);
    const real C_TA = std::pow(10.0, -std::floor(beta_A));
    
    // 计算δ_k
    const int delta0 = (chi0 < C_TA) ? 0 : 1;
    const int delta1 = (chi1 < C_TA) ? 0 : 1;
    const int delta2 = (chi2 < C_TA) ? 0 : 1;

    // 根据δ值选择系数
    real a_im2, a_im1, a_i, a_ip1, a_ip2;
    
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a_im2 = 3.0/128.0;
        a_im1 = -5.0/32.0;
        a_i = 45.0/64.0;
        a_ip1 = 15.0/32.0;
        a_ip2 = -5.0/128.0;
    }
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a_im2 = 0.0;
        a_im1 = -1.0/16.0;
        a_i = 9.0/16.0;
        a_ip1 = 9.0/16.0;
        a_ip2 = -1.0/16.0;
    }
    else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
             (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a_im2 = 0.0;
        a_im1 = 0.0;
        a_i = 3.0/8.0;
        a_ip1 = 3.0/4.0;
        a_ip2 = -1.0/8.0;
    }
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a_im2 = 1.0/16.0;
        a_im1 = -5.0/16.0;
        a_i = 15.0/16.0;
        a_ip1 = 5.0/16.0;
        a_ip2 = 0.0;
    }
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a_im2 = 0.0;
        a_im1 = -1.0/8.0;
        a_i = 3.0/4.0;
        a_ip1 = 3.0/8.0;
        a_ip2 = 0.0;
    }
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a_im2 = 3.0/8.0;
        a_im1 = -5.0/4.0;
        a_i = 15.0/8.0;
        a_ip1 = 0.0;
        a_ip2 = 0.0;
    }
    else {  // delta0 == 0 && delta1 == 0 && delta2 == 0
        a_im2 = 0.0;
        a_im1 = 0.0;
        a_i = 1.0;
        a_ip1 = 0.0;
        a_ip2 = 0.0;
    }

    real mmmmmmmmm1 = 0.0;
    real mmmmmmmmm2 = 0.0;

    // 计算最终插值
    return a_im2 * f_im2 + a_im1 * f_im1 + a_i * f_i + a_ip1 * f_ip1 + a_ip2 * f_ip2;
}

inline real whf_ai_TcnsN_A_2(std::array<real, 5> q) {
    // 常量定义
    constexpr real C = 1.0;
    constexpr real q_exp = 6.0;  // 指数q
    constexpr real epsilon = 1e-40;
    constexpr real Cr = 0.24;
    constexpr real xi = 1e-3;
    constexpr real epsilon_A = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;
    constexpr real alpha1 = 10.0;  // 示例值，需根据实际调整
    constexpr real alpha2 = 5.0;  // 示例值，需根据实际调整
    
    // 提取输入值
    real q_im2 = q[0];  // q_{i-2}
    real q_im1 = q[1];  // q_{i-1}
    real q_i   = q[2];  // q_i
    real q_ip1 = q[3];  // q_{i+1}
    real q_ip2 = q[4];  // q_{i+2}
    
    // 计算光滑度量 β_k
    real beta0 = std::pow(q_im2 - 2*q_im1 + q_i, 2) 
               + 0.25 * std::pow(q_im2 - 4*q_im1 + 3*q_i, 2);
               
    real beta1 = std::pow(q_im1 - 2*q_i + q_ip1, 2) 
               + 0.25 * std::pow(q_im1 - q_ip1, 2);
               
    real beta2 = std::pow(q_i - 2*q_ip1 + q_ip2, 2) 
               + 0.25 * std::pow(3*q_i - 4*q_ip1 + q_ip2, 2);
    
    // 计算全局光滑因子 τ
    real tau = std::abs(beta2 - beta0);
    
    // 计算 γ_k
    real gamma0 = std::pow(C + tau/(beta0 + epsilon), q_exp);
    real gamma1 = std::pow(C + tau/(beta1 + epsilon), q_exp);
    real gamma2 = std::pow(C + tau/(beta2 + epsilon), q_exp);
    
    // 归一化光滑度量 χ_k
    real sum_gamma = gamma0 + gamma1 + gamma2;
    real chi0 = gamma0 / sum_gamma;
    real chi1 = gamma1 / sum_gamma;
    real chi2 = gamma2 / sum_gamma;
    
    // 计算 Δq
    real Dq_im32 = q_im1 - q_im2;  // Δq_{i-3/2}
    real Dq_im12 = q_i - q_im1;   // Δq_{i-1/2}
    real Dq_ip12 = q_ip1 - q_i;   // Δq_{i+1/2}
    real Dq_ip32 = q_ip2 - q_ip1; // Δq_{i+3/2}
    
    // 计算 η
    real eta_im1 = (std::abs(2*Dq_im12*Dq_im32) + epsilon_A) 
                 / (Dq_im12*Dq_im12 + Dq_im32*Dq_im32 + epsilon_A);
                 
    real eta_i = (std::abs(2*Dq_ip12*Dq_im12) + epsilon_A) 
               / (Dq_ip12*Dq_ip12 + Dq_im12*Dq_im12 + epsilon_A);
               
    real eta_ip1 = (std::abs(2*Dq_ip32*Dq_ip12) + epsilon_A) 
                 / (Dq_ip32*Dq_ip32 + Dq_ip12*Dq_ip12 + epsilon_A);
    
    // 计算 η_{i+1/2}
    real eta_half = std::min({eta_im1, eta_i, eta_ip1});
    
    // 计算 m
    real m = 1.0 - std::min(1.0, eta_half/Cr);
    
    // 计算 g(m)
    real g_m = std::pow(1 - m, 4) * (1 + 4*m);
    
    // 计算 β_A 和 C_TA
    real beta_A = alpha1 - alpha2*(1 - g_m);
    real C_TA = std::pow(10, -std::floor(beta_A));
    
    // 计算 δ_k
    int delta0 = (chi0 < C_TA) ? 0 : 1;
    int delta1 = (chi1 < C_TA) ? 0 : 1;
    int delta2 = (chi2 < C_TA) ? 0 : 1;
    
    // 根据δ_k选择插值系数
    real a_im2, a_im1, a_i, a_ip1, a_ip2;
    
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a_im2 = 3.0/128; a_im1 = -5.0/32;
        a_i = 45.0/64; a_ip1 = 15.0/32; a_ip2 = -5.0/128;
    } 
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a_im2 = 0; a_im1 = -1.0/16;
        a_i = 9.0/16; a_ip1 = 9.0/16; a_ip2 = -1.0/16;
    } 
    else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) ||
             (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a_im2 = 0; a_im1 = 0;
        a_i = 3.0/8; a_ip1 = 3.0/4; a_ip2 = -1.0/8;
    } 
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a_im2 = 1.0/16; a_im1 = -5.0/16;
        a_i = 15.0/16; a_ip1 = 5.0/16; a_ip2 = 0;
    } 
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a_im2 = 0; a_im1 = -1.0/8;
        a_i = 3.0/4; a_ip1 = 3.0/8; a_ip2 = 0;
    } 
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a_im2 = 3.0/8; a_im1 = -5.0/4;
        a_i = 15.0/8; a_ip1 = 0; a_ip2 = 0;
    } 
    else {  // (0,0,0)
        a_im2 = 0; a_im1 = 0;
        a_i = 1; a_ip1 = 0; a_ip2 = 0;
    }
    
    // 计算最终插值结果
    return a_im2 * q_im2 + a_im1 * q_im1 + a_i * q_i + a_ip1 * q_ip1 + a_ip2 * q_ip2;
}

inline real whf_ai_TcnsN_AS_1(std::array<real, 5> q) {
    // 常量定义
    constexpr real xi = 1e-3;
    constexpr real Cr = 0.24;
    constexpr real alpha1 = 10.0;
    constexpr real alpha2 = 5.0;
    constexpr real qA = 6.0;
    constexpr real C = 1.0;
    constexpr real epsilon_A = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;

    // 解包输入值
    real q_im2 = q[0]; // q_{i-2}
    real q_im1 = q[1]; // q_{i-1}
    real q_i   = q[2]; // q_i
    real q_ip1 = q[3]; // q_{i+1}
    real q_ip2 = q[4]; // q_{i+2}

    // 计算差分值
    real dq_im3_2 = q_im1 - q_im2; // Δq_{i-3/2}
    real dq_im1_2 = q_i - q_im1;   // Δq_{i-1/2}
    real dq_ip1_2 = q_ip1 - q_i;   // Δq_{i+1/2}
    real dq_ip3_2 = q_ip2 - q_ip1; // Δq_{i+3/2}

    // 计算局部光滑因子
    real beta0 = std::pow(q_im2 - 2*q_im1 + q_i, 2) + 
                 0.25 * std::pow(q_im2 - 4*q_im1 + 3*q_i, 2);
    
    real beta1 = std::pow(q_im1 - 2*q_i + q_ip1, 2) + 
                 0.25 * std::pow(q_im1 - q_ip1, 2);
    
    real beta2 = std::pow(q_i - 2*q_ip1 + q_ip2, 2) + 
                 0.25 * std::pow(3*q_i - 4*q_ip1 + q_ip2, 2);

    // 计算全局光滑因子
    real tau = std::abs(beta2 - beta0);

    // 计算η值
    real eta_im1 = (std::abs(2*dq_im1_2*dq_im3_2) + epsilon_A) / 
                   (std::pow(dq_im1_2, 2) + std::pow(dq_im3_2, 2) + epsilon_A);
    
    real eta_i = (std::abs(2*dq_ip1_2*dq_im1_2) + epsilon_A) / 
                 (std::pow(dq_ip1_2, 2) + std::pow(dq_im1_2, 2) + epsilon_A);
    
    real eta_ip1 = (std::abs(2*dq_ip3_2*dq_ip1_2) + epsilon_A) / 
                   (std::pow(dq_ip3_2, 2) + std::pow(dq_ip1_2, 2) + epsilon_A);
    
    real eta_min = std::min({eta_im1, eta_i, eta_ip1});

    // 计算m和g(m)
    real m = 1.0 - std::min(1.0, eta_min / Cr);
    real g_m = std::pow(1.0 - m, 4) * (1.0 + 4.0*m);
    
    // 计算自适应参数
    real beta_A = alpha1 - alpha2 * (1.0 - g_m);
    real C_TA = std::pow(10.0, -std::floor(beta_A));
    real C_T_prime = std::pow(1.5 * C_TA / (1.0 - C_TA), 1.0/qA);

    // 确定β_j
    real beta_j = std::min({beta0, beta1, beta2});

    // 计算δ_k
    int delta0 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta0) ? 0 : 1;
    int delta1 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta1) ? 0 : 1;
    int delta2 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta2) ? 0 : 1;

    // 根据δ值选择插值系数
    real a0, a1, a2, a3, a4;
    
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a0 = 3.0/128; a1 = -5.0/32; a2 = 45.0/64; a3 = 15.0/32; a4 = -5.0/128;
    } 
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a0 = 0; a1 = -1.0/16; a2 = 9.0/16; a3 = 9.0/16; a4 = -1.0/16;
    } 
    else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
             (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a0 = 0; a1 = 0; a2 = 3.0/8; a3 = 3.0/4; a4 = -1.0/8;
    } 
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a0 = 1.0/16; a1 = -5.0/16; a2 = 15.0/16; a3 = 5.0/16; a4 = 0;
    } 
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a0 = 0; a1 = -1.0/8; a2 = 3.0/4; a3 = 3.0/8; a4 = 0;
    } 
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a0 = 3.0/8; a1 = -5.0/4; a2 = 15.0/8; a3 = 0; a4 = 0;
    } 
    else { // (0,0,0)
        a0 = 0; a1 = 0; a2 = 1; a3 = 0; a4 = 0;
    }

    // 计算最终插值结果
    return a0*q_im2 + a1*q_im1 + a2*q_i + a3*q_ip1 + a4*q_ip2;
}

inline real whf_ai_TcnsN_AS_2(std::array<real, 5> q) {
    // 提取输入值
    real q_im2 = q[0]; // q_{i-2}
    real q_im1 = q[1]; // q_{i-1}
    real q_i   = q[2]; // q_i
    real q_ip1 = q[3]; // q_{i+1}
    real q_ip2 = q[4]; // q_{i+2}

    // 常数定义
    const real xi = 1e-3;      // ξ
    const real Cr = 0.24;      // C_r
    const real alpha1 = 10.0;  // α1
    const real alpha2 = 5.0;   // α2
    const real qA = 6.0;       // q_A

    // 计算差分值
    real dq_im3_2 = q_im1 - q_im2; // Δq_{i-3/2}
    real dq_im1_2 = q_i - q_im1;   // Δq_{i-1/2}
    real dq_ip1_2 = q_ip1 - q_i;   // Δq_{i+1/2}
    real dq_ip3_2 = q_ip2 - q_ip1; // Δq_{i+3/2}

    // 计算局部光滑因子β
    real beta0 = std::pow(q_im2 - 2*q_im1 + q_i, 2) + 
                 0.25 * std::pow(q_im2 - 4*q_im1 + 3*q_i, 2);
    real beta1 = std::pow(q_im1 - 2*q_i + q_ip1, 2) + 
                 0.25 * std::pow(q_im1 - q_ip1, 2);
    real beta2 = std::pow(q_i - 2*q_ip1 + q_ip2, 2) + 
                 0.25 * std::pow(3*q_i - 4*q_ip1 + q_ip2, 2);

    // 计算全局光滑因子τ
    real tau = std::abs(beta2 - beta0);

    // 计算η值
    real epsA = (0.9 * Cr / (1 - 0.9 * Cr)) * xi * xi;
    
    real eta_im1 = (std::abs(2 * dq_im1_2 * dq_im3_2) + epsA) /
                   (dq_im1_2*dq_im1_2 + dq_im3_2*dq_im3_2 + epsA);
    real eta_i = (std::abs(2 * dq_ip1_2 * dq_im1_2) + epsA) /
                 (dq_ip1_2*dq_ip1_2 + dq_im1_2*dq_im1_2 + epsA);
    real eta_ip1 = (std::abs(2 * dq_ip3_2 * dq_ip1_2) + epsA) /
                   (dq_ip3_2*dq_ip3_2 + dq_ip1_2*dq_ip1_2 + epsA);
    
    real eta_ip1_2 = std::min({eta_im1, eta_i, eta_ip1});

    // 计算m和g(m)
    real m = 1.0 - std::min(1.0, eta_ip1_2 / Cr);
    real gm = std::pow(1.0 - m, 4) * (1.0 + 4.0 * m);

    // 计算βA和C_TA
    real betaA = alpha1 - alpha2 * (1.0 - gm);
    real C_TA = std::pow(10.0, -std::floor(betaA));

    // 计算辅助截断函数C_T'
    real C_T_prime = std::pow(1.5 * C_TA / (1.0 - C_TA), 1.0 / qA);

    // 计算βj = min{β0, β1, β2}
    real beta_j = std::min({beta0, beta1, beta2});

    // 计算δk (k=0,1,2)
    real C = 1.0;
    real delta0 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta0) ? 0 : 1;
    real delta1 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta1) ? 0 : 1;
    real delta2 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta2) ? 0 : 1;

    // 根据δ值选择插值系数
    real a_im2, a_im1, a_i, a_ip1, a_ip2;

    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a_im2 = 3.0/128.0; a_im1 = -5.0/32.0; a_i = 45.0/64.0; a_ip1 = 15.0/32.0; a_ip2 = -5.0/128.0;
    } else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a_im2 = 0.0; a_im1 = -1.0/16.0; a_i = 9.0/16.0; a_ip1 = 9.0/16.0; a_ip2 = -1.0/16.0;
    } else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
               (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a_im2 = 0.0; a_im1 = 0.0; a_i = 3.0/8.0; a_ip1 = 3.0/4.0; a_ip2 = -1.0/8.0;
    } else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a_im2 = 1.0/16.0; a_im1 = -5.0/16.0; a_i = 15.0/16.0; a_ip1 = 5.0/16.0; a_ip2 = 0.0;
    } else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a_im2 = 0.0; a_im1 = -1.0/8.0; a_i = 3.0/4.0; a_ip1 = 3.0/8.0; a_ip2 = 0.0;
    } else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a_im2 = 3.0/8.0; a_im1 = -5.0/4.0; a_i = 15.0/8.0; a_ip1 = 0.0; a_ip2 = 0.0;
    } else { // delta0==0 && delta1==0 && delta2==0
        a_im2 = 0.0; a_im1 = 0.0; a_i = 1.0; a_ip1 = 0.0; a_ip2 = 0.0;
    }

    // 计算最终插值结果
    return a_im2 * q_im2 + a_im1 * q_im1 + a_i * q_i + a_ip1 * q_ip1 + a_ip2 * q_ip2;
}

inline real whf_ai_TcnsN_LADS_yuanbao(std::array<real, 5> q) {
    // 常量定义
    const real epsilon = 1e-40;
    const real H = 10.0;
    const int B_l = 5;
    const int B_u = 10;
    const int q_A = 6;
    const real C = 1.0;
    
    // 提取输入值
    real q_i2 = q[0]; // q_{i-2}
    real q_i1 = q[1]; // q_{i-1}
    real q_i = q[2];  // q_i
    real q_i_1 = q[3]; // q_{i+1}
    real q_i_2 = q[4]; // q_{i+2}
    
    // 计算局部光滑因子β
    real beta0 = std::pow(q_i2 - 2*q_i1 + q_i, 2) + 
                0.25 * std::pow(q_i2 - 4*q_i1 + 3*q_i, 2);
    
    real beta1 = std::pow(q_i1 - 2*q_i + q_i_1, 2) + 
                0.25 * std::pow(q_i1 - q_i_1, 2);
    
    real beta2 = std::pow(q_i - 2*q_i_1 + q_i_2, 2) + 
                0.25 * std::pow(3*q_i - 4*q_i_1 + q_i_2, 2);
    
    // 计算全局光滑因子τ
    real tau = std::abs(beta2 - beta0);
    
    // 计算χ_k = τ/(β_k + ε)
    real chi0 = tau / (beta0 + epsilon);
    real chi1 = tau / (beta1 + epsilon);
    real chi2 = tau / (beta2 + epsilon);
    
    // 找到最大的χ_k
    real max_chi = std::max({chi0, chi1, chi2});
    
    // 计算θ
    real theta = 1.0 / (1.0 + max_chi / H);
    
    // 计算m和C_T
    int m = B_l + static_cast<int>(std::floor(theta * (B_u - B_l)));
    real C_T = std::pow(10.0, -m);
    
    // 计算辅助截断函数C_T'
    real C_T_prime = std::pow(1.5 * C_T / (1.0 - C_T), 1.0 / q_A);
    
    // 找到最小的β_j
    real beta_j = std::min({beta0, beta1, beta2});
    
    // 计算δ_k
    int delta0, delta1, delta2;
    
    // 计算条件表达式
    real condition_expr = C_T_prime * tau - C * (1.0 - C_T_prime) * beta_j;
    
    // 计算δ0
    if (beta_j * tau < condition_expr * beta0) {
        delta0 = 0;
    } else {
        delta0 = 1;
    }
    
    // 计算δ1
    if (beta_j * tau < condition_expr * beta1) {
        delta1 = 0;
    } else {
        delta1 = 1;
    }
    
    // 计算δ2
    if (beta_j * tau < condition_expr * beta2) {
        delta2 = 0;
    } else {
        delta2 = 1;
    }
    
    // 根据δ值选择插值系数
    real a_i2, a_i1, a_i, a_i1_, a_i2_; // a_{i-2}, a_{i-1}, a_i, a_{i+1}, a_{i+2}
    
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a_i2 = 3.0/128.0;
        a_i1 = -5.0/32.0;
        a_i = 45.0/64.0;
        a_i1_ = 15.0/32.0;
        a_i2_ = -5.0/128.0;
    } 
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a_i2 = 0.0;
        a_i1 = -1.0/16.0;
        a_i = 9.0/16.0;
        a_i1_ = 9.0/16.0;
        a_i2_ = -1.0/16.0;
    } 
    else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
             (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a_i2 = 0.0;
        a_i1 = 0.0;
        a_i = 3.0/8.0;
        a_i1_ = 3.0/4.0;
        a_i2_ = -1.0/8.0;
    } 
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a_i2 = 1.0/16.0;
        a_i1 = -5.0/16.0;
        a_i = 15.0/16.0;
        a_i1_ = 5.0/16.0;
        a_i2_ = 0.0;
    } 
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a_i2 = 0.0;
        a_i1 = -1.0/8.0;
        a_i = 3.0/4.0;
        a_i1_ = 3.0/8.0;
        a_i2_ = 0.0;
    } 
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a_i2 = 3.0/8.0;
        a_i1 = -5.0/4.0;
        a_i = 15.0/8.0;
        a_i1_ = 0.0;
        a_i2_ = 0.0;
    } 
    else { // (0,0,0)
        a_i2 = 0.0;
        a_i1 = 0.0;
        a_i = 1.0;
        a_i1_ = 0.0;
        a_i2_ = 0.0;
    }
    
    // 计算最终的插值
    return a_i2 * q_i2 + a_i1 * q_i1 + a_i * q_i + a_i1_ * q_i_1 + a_i2_ * q_i_2;
}


inline real whf_ai_TcnsN_LADS_doubao(std::array<real, 5> q) {
    // 提取输入数据
    real q_im2 = q[0];  // q_{i-2}
    real q_im1 = q[1];  // q_{i-1}
    real q_i   = q[2];  // q_i
    real q_ip1 = q[3];  // q_{i+1}
    real q_ip2 = q[4];  // q_{i+2}

    // 计算局部光滑因子 beta0, beta1, beta2
    real beta0 = std::pow(q_im2 - 2*q_im1 + q_i, 2) + 0.25 * std::pow(q_im2 - 4*q_im1 + 3*q_i, 2);
    real beta1 = std::pow(q_im1 - 2*q_i + q_ip1, 2) + 0.25 * std::pow(q_im1 - q_ip1, 2);
    real beta2 = std::pow(q_i - 2*q_ip1 + q_ip2, 2) + 0.25 * std::pow(3*q_i - 4*q_ip1 + q_ip2, 2);

    // 计算全局光滑因子 tau
    real tau = std::abs(beta2 - beta0);

    // 计算辅助参数
    const real eps = 1e-40;
    real chi0 = tau / (beta0 + eps);
    real chi1 = tau / (beta1 + eps);
    real chi2 = tau / (beta2 + eps);
    real max_chi = std::max({chi0, chi1, chi2});

    const real H = 10.0;
    real theta = 1.0 / (1.0 + (max_chi / H));

    const int Bu = 10;
    const int Bl = 5;
    int m = Bl + static_cast<int>(std::floor(theta * (Bu - Bl)));
    real CT = std::pow(10.0, -m);

    // 计算辅助截断函数 C_T'
    const int qA = 6;
    real CT_prime = std::pow((1.5 * CT) / (1.0 - CT), 1.0 / qA);

    // 确定 beta_j (最小的 beta)
    real betaj = std::min({beta0, beta1, beta2});
    const real C = 1.0;

    // 计算截断函数 delta0, delta1, delta2
    auto compute_delta = [&](real betak) {
        real rhs = (CT_prime * tau - C * (1.0 - CT_prime) * betaj) * betak;
        return (betaj * tau < rhs) ? 0 : 1;
    };

    int delta0 = compute_delta(beta0);
    int delta1 = compute_delta(beta1);
    int delta2 = compute_delta(beta2);

    // 根据 delta 组合确定系数 a
    real a_im2, a_im1, a_i, a_ip1, a_ip2;

    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a_im2 = 3.0 / 128.0;
        a_im1 = -5.0 / 32.0;
        a_i   = 45.0 / 64.0;
        a_ip1 = 15.0 / 32.0;
        a_ip2 = -5.0 / 128.0;
    } else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a_im2 = 0.0;
        a_im1 = -1.0 / 16.0;
        a_i   = 9.0 / 16.0;
        a_ip1 = 9.0 / 16.0;
        a_ip2 = -1.0 / 16.0;
    } else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
               (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a_im2 = 0.0;
        a_im1 = 0.0;
        a_i   = 3.0 / 8.0;
        a_ip1 = 3.0 / 4.0;
        a_ip2 = -1.0 / 8.0;
    } else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a_im2 = 1.0 / 16.0;
        a_im1 = -5.0 / 16.0;
        a_i   = 15.0 / 16.0;
        a_ip1 = 5.0 / 16.0;
        a_ip2 = 0.0;
    } else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a_im2 = 0.0;
        a_im1 = -1.0 / 8.0;
        a_i   = 3.0 / 4.0;
        a_ip1 = 3.0 / 8.0;
        a_ip2 = 0.0;
    } else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a_im2 = 3.0 / 8.0;
        a_im1 = -5.0 / 4.0;
        a_i   = 15.0 / 8.0;
        a_ip1 = 0.0;
        a_ip2 = 0.0;
    } else {  // delta0 == 0 && delta1 == 0 && delta2 == 0
        a_im2 = 0.0;
        a_im1 = 0.0;
        a_i   = 1.0;
        a_ip1 = 0.0;
        a_ip2 = 0.0;
    }

    // 计算并返回最终插值结果
    return a_im2 * q_im2 + a_im1 * q_im1 + a_i * q_i + a_ip1 * q_ip1 + a_ip2 * q_ip2;
}

inline real whf_ai_TcnsN_AZS_doubao(std::array<real, 5> q) {
    // 提取输入数组中的各分量，对应q_{i-2}到q_{i+2}
    real q0 = q[0];  // q_{i-2}
    real q1 = q[1];  // q_{i-1}
    real q2 = q[2];  // q_i
    real q3 = q[3];  // q_{i+1}
    real q4 = q[4];  // q_{i+2}

    // 计算局部光滑因子β0、β1、β2
    real term1 = q0 - 2 * q1 + q2;
    real term2 = q0 - 4 * q1 + 3 * q2;
    real beta0 = term1 * term1 + 0.25 * term2 * term2;

    term1 = q1 - 2 * q2 + q3;
    term2 = q1 - q3;
    real beta1 = term1 * term1 + 0.25 * term2 * term2;

    term1 = q2 - 2 * q3 + q4;
    term2 = 3 * q2 - 4 * q3 + q4;
    real beta2 = term1 * term1 + 0.25 * term2 * term2;

    // 计算全局光滑因子τ
    real tau = std::abs(beta2 - beta0);

    // 计算C_T相关参数
    real tau_A = std::abs(beta0 + beta2 - 2 * beta1);
    const real lambda_d = 1e-5;
    real C_T = std::min(lambda_d, tau_A);

    // 计算辅助截断函数C_T'
    const real q_A = 6.0;
    real numerator = 1.5 * C_T;
    real denominator = 1.0 - C_T;
    real ratio = numerator / denominator;
    real C_T_prime = std::pow(ratio, 1.0 / q_A);

    // 确定β_j（β0、β1、β2中的最小值）
    real beta_j = std::min({beta0, beta1, beta2});

    // 计算δ0、δ1、δ2（使用lambda函数简化重复逻辑）
    auto compute_delta = [&](real beta_k) {
        real left = beta_j * tau;
        real right = (C_T_prime * tau - (1.0 - C_T_prime) * beta_j) * beta_k;
        return (left < right) ? 0 : 1;
    };

    int delta0 = compute_delta(beta0);
    int delta1 = compute_delta(beta1);
    int delta2 = compute_delta(beta2);

    // 根据δ组合确定插值系数a_{i-2}到a_{i+2}
    real a0, a1, a2_coeff, a3, a4;  // 用a2_coeff避免与q2重名
    int delta = (delta0 << 2) | (delta1 << 1) | delta2;

    switch (delta) {
        case 7:  // δ0=1, δ1=1, δ2=1
            a0 = 3.0 / 128.0;
            a1 = -5.0 / 32.0;
            a2_coeff = 45.0 / 64.0;
            a3 = 15.0 / 32.0;
            a4 = -5.0 / 128.0;
            break;
        case 3:  // δ0=0, δ1=1, δ2=1
            a0 = 0.0;
            a1 = -1.0 / 16.0;
            a2_coeff = 9.0 / 16.0;
            a3 = 9.0 / 16.0;
            a4 = -1.0 / 16.0;
            break;
        case 5:  // δ0=1, δ1=0, δ2=1
        case 1:  // δ0=0, δ1=0, δ2=1
            a0 = 0.0;
            a1 = 0.0;
            a2_coeff = 3.0 / 8.0;
            a3 = 3.0 / 4.0;
            a4 = -1.0 / 8.0;
            break;
        case 6:  // δ0=1, δ1=1, δ2=0
            a0 = 1.0 / 16.0;
            a1 = -5.0 / 16.0;
            a2_coeff = 15.0 / 16.0;
            a3 = 5.0 / 16.0;
            a4 = 0.0;
            break;
        case 2:  // δ0=0, δ1=1, δ2=0
            a0 = 0.0;
            a1 = -1.0 / 8.0;
            a2_coeff = 3.0 / 4.0;
            a3 = 3.0 / 8.0;
            a4 = 0.0;
            break;
        case 4:  // δ0=1, δ1=0, δ2=0
            a0 = 3.0 / 8.0;
            a1 = -5.0 / 4.0;
            a2_coeff = 15.0 / 8.0;
            a3 = 0.0;
            a4 = 0.0;
            break;
        case 0:  // δ0=0, δ1=0, δ2=0
        default:
            a0 = 0.0;
            a1 = 0.0;
            a2_coeff = 1.0;
            a3 = 0.0;
            a4 = 0.0;
            break;
    }

    // 计算并返回最终结果\hat{f}_{i+1/2}
    return a0 * q0 + a1 * q1 + a2_coeff * q2 + a3 * q3 + a4 * q4;
}


inline real whf_ai_TcnsN_AZS_yuanbao(std::array<real, 5> q) {
    // 解包输入值
    real q_i2 = q[0]; // q_{i-2}
    real q_i1 = q[1]; // q_{i-1}
    real q_i  = q[2]; // q_i
    real q_i1p = q[3]; // q_{i+1}
    real q_i2p = q[4]; // q_{i+2}

    // 计算局部光滑因子β
    real beta0 = std::pow(q_i2 - 2*q_i1 + q_i, 2) + 
                 0.25 * std::pow(q_i2 - 4*q_i1 + 3*q_i, 2);
    
    real beta1 = std::pow(q_i1 - 2*q_i + q_i1p, 2) + 
                 0.25 * std::pow(q_i1 - q_i1p, 2);
    
    real beta2 = std::pow(q_i - 2*q_i1p + q_i2p, 2) + 
                 0.25 * std::pow(3*q_i - 4*q_i1p + q_i2p, 2);

    // 计算全局光滑因子τ
    real tau = std::fabs(beta2 - beta0);

    // 计算自适应参数
    const real lambda_d = 1e-5;
    real tau_A = std::fabs(beta0 + beta2 - 2*beta1);
    real C_T = std::min(lambda_d, tau_A);

    // 计算辅助截断函数C_T'
    const real q_A = 6.0;
    real C_T_prime = std::pow(1.5 * C_T / (1 - C_T), 1.0/q_A);

    // 计算最小光滑因子β_j
    real beta_j = std::min({beta0, beta1, beta2});

    // 计算截断函数δ
    const real C = 1.0;
    real delta0 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta0) ? 0 : 1;
    real delta1 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta1) ? 0 : 1;
    real delta2 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta2) ? 0 : 1;

    // 定义系数表 [a_{i-2}, a_{i-1}, a_i, a_{i+1}, a_{i+2}]
    const std::array<std::array<real, 5>, 8> coefficients = {{
        { 3.0/128, -5.0/32, 45.0/64, 15.0/32, -5.0/128 }, // δ=(1,1,1)
        { 0.0,     -1.0/16,  9.0/16,  9.0/16, -1.0/16  }, // δ=(0,1,1)
        { 0.0,      0.0,      3.0/8,   3.0/4,  -1.0/8  }, // δ=(1,0,1)
        { 0.0,      0.0,      3.0/8,   3.0/4,  -1.0/8  }, // δ=(0,0,1)
        { 1.0/16,  -5.0/16, 15.0/16,  5.0/16,  0.0     }, // δ=(1,1,0)
        { 0.0,     -1.0/8,    3.0/4,   3.0/8,   0.0     }, // δ=(0,1,0)
        { 3.0/8,   -5.0/4,   15.0/8,   0.0,      0.0     }, // δ=(1,0,0)
        { 0.0,      0.0,      1.0,     0.0,      0.0     }  // δ=(0,0,0)
    }};

    // 根据δ值选择系数
    size_t index = 0;
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) index = 0;
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) index = 1;
    else if (delta0 == 1 && delta1 == 0 && delta2 == 1) index = 2;
    else if (delta0 == 0 && delta1 == 0 && delta2 == 1) index = 3;
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) index = 4;
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) index = 5;
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) index = 6;
    else index = 7; // (0,0,0)

    // 计算插值结果
    const auto& a = coefficients[index];
    return a[0]*q_i2 + a[1]*q_i1 + a[2]*q_i + a[3]*q_i1p + a[4]*q_i2p;
}




//【王鸿飞】end插值格式开发（2.0）AI-原格式

//【王鸿飞】begin插值格式开发（3.0）AI-新格式

inline real whf_ai_TcnsN_myASF002_1(std::array<real, 5> q) {
    // 常量定义
    constexpr real xi = 1e-3;
    constexpr real Cr = 0.24;
    // constexpr real alpha1 = 11.0;
    // constexpr real alpha2 = 6.0;
    // constexpr real qA = 6.0;
    constexpr real C = 1.0;
    constexpr real epsilon_A = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;

    // 解包输入值
    real q_im2 = q[0]; // q_{i-2}
    real q_im1 = q[1]; // q_{i-1}
    real q_i   = q[2]; // q_i
    real q_ip1 = q[3]; // q_{i+1}
    real q_ip2 = q[4]; // q_{i+2}

    // 计算差分值
    real dq_im3_2 = q_im1 - q_im2; // Δq_{i-3/2}
    real dq_im1_2 = q_i - q_im1;   // Δq_{i-1/2}
    real dq_ip1_2 = q_ip1 - q_i;   // Δq_{i+1/2}
    real dq_ip3_2 = q_ip2 - q_ip1; // Δq_{i+3/2}

    // 计算局部光滑因子
    real beta0 = std::pow(q_im2 - 2*q_im1 + q_i, 2) + 
                 0.25 * std::pow(q_im2 - 4*q_im1 + 3*q_i, 2);
    
    real beta1 = std::pow(q_im1 - 2*q_i + q_ip1, 2) + 
                 0.25 * std::pow(q_im1 - q_ip1, 2);
    
    real beta2 = std::pow(q_i - 2*q_ip1 + q_ip2, 2) + 
                 0.25 * std::pow(3*q_i - 4*q_ip1 + q_ip2, 2);

    // 计算全局光滑因子
    real tau = std::abs(beta2 - beta0);

    // 计算η值
    real eta_im1 = (std::abs(2.0*dq_im1_2*dq_im3_2) + epsilon_A) / 
                   (std::pow(dq_im1_2, 2) + std::pow(dq_im3_2, 2) + epsilon_A);
    
    real eta_i = (std::abs(2.0*dq_ip1_2*dq_im1_2) + epsilon_A) / 
                 (std::pow(dq_ip1_2, 2) + std::pow(dq_im1_2, 2) + epsilon_A);
    
    real eta_ip1 = (std::abs(2.0*dq_ip3_2*dq_ip1_2) + epsilon_A) / 
                   (std::pow(dq_ip3_2, 2) + std::pow(dq_ip1_2, 2) + epsilon_A);
    
    real eta_min = std::min({eta_im1, eta_i, eta_ip1});

    // 计算min
    real min = std::min(1.0, eta_min / Cr);
    
    // 计算C_T_prime
    real C_T_prime = 0.157/(11.351092*min*min - 5.740650*min + 1.506631);

    // 确定β_j
    real beta_j = std::min({beta0, beta1, beta2});

    // 计算δ_k
    int delta0 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta0) ? 0 : 1;
    int delta1 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta1) ? 0 : 1;
    int delta2 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta2) ? 0 : 1;

    // 根据δ值选择插值系数
    real a0, a1, a2, a3, a4;
    
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a0 = 3.0/128; a1 = -5.0/32; a2 = 45.0/64; a3 = 15.0/32; a4 = -5.0/128;
    } 
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a0 = 0; a1 = -1.0/16; a2 = 9.0/16; a3 = 9.0/16; a4 = -1.0/16;
    } 
    else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
             (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a0 = 0; a1 = 0; a2 = 3.0/8; a3 = 3.0/4; a4 = -1.0/8;
    } 
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a0 = 1.0/16; a1 = -5.0/16; a2 = 15.0/16; a3 = 5.0/16; a4 = 0;
    } 
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a0 = 0; a1 = -1.0/8; a2 = 3.0/4; a3 = 3.0/8; a4 = 0;
    } 
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a0 = 3.0/8; a1 = -5.0/4; a2 = 15.0/8; a3 = 0; a4 = 0;
    } 
    else { // (0,0,0)
        a0 = 0; a1 = 0; a2 = 1; a3 = 0; a4 = 0;
    }

    // 计算最终插值结果
    return a0*q_im2 + a1*q_im1 + a2*q_i + a3*q_ip1 + a4*q_ip2;
}


inline real whf_ai_TcnsN_myASF002_2(std::array<real, 5> q) {
    // 提取输入值
    real q_im2 = q[0]; // q_{i-2}
    real q_im1 = q[1]; // q_{i-1}
    real q_i   = q[2]; // q_i
    real q_ip1 = q[3]; // q_{i+1}
    real q_ip2 = q[4]; // q_{i+2}

    // 常数定义
    const real xi = 1e-3;      // ξ
    const real Cr = 0.24;      // C_r
    // const real alpha1 = 10.0;  // α1
    // const real alpha2 = 5.0;   // α2
    // const real qA = 6.0;       // q_A

    // 计算差分值
    real dq_im3_2 = q_im1 - q_im2; // Δq_{i-3/2}
    real dq_im1_2 = q_i - q_im1;   // Δq_{i-1/2}
    real dq_ip1_2 = q_ip1 - q_i;   // Δq_{i+1/2}
    real dq_ip3_2 = q_ip2 - q_ip1; // Δq_{i+3/2}

    // 计算局部光滑因子β
    real beta0 = std::pow(q_im2 - 2*q_im1 + q_i, 2) + 
                 0.25 * std::pow(q_im2 - 4*q_im1 + 3*q_i, 2);
    real beta1 = std::pow(q_im1 - 2*q_i + q_ip1, 2) + 
                 0.25 * std::pow(q_im1 - q_ip1, 2);
    real beta2 = std::pow(q_i - 2*q_ip1 + q_ip2, 2) + 
                 0.25 * std::pow(3*q_i - 4*q_ip1 + q_ip2, 2);

    // 计算全局光滑因子τ
    real tau = std::abs(beta2 - beta0);

    // 计算η值
    real epsA = (0.9 * Cr / (1 - 0.9 * Cr)) * xi * xi;
    
    real eta_im1 = (std::abs(2 * dq_im1_2 * dq_im3_2) + epsA) /
                   (dq_im1_2*dq_im1_2 + dq_im3_2*dq_im3_2 + epsA);
    real eta_i = (std::abs(2 * dq_ip1_2 * dq_im1_2) + epsA) /
                 (dq_ip1_2*dq_ip1_2 + dq_im1_2*dq_im1_2 + epsA);
    real eta_ip1 = (std::abs(2 * dq_ip3_2 * dq_ip1_2) + epsA) /
                   (dq_ip3_2*dq_ip3_2 + dq_ip1_2*dq_ip1_2 + epsA);
    
    real eta_ip1_2 = std::min({eta_im1, eta_i, eta_ip1});

    // 计算min
    real min = std::min(1.0, eta_ip1_2 / Cr);

    // 计算辅助截断函数C_T'
    real C_T_prime = 0.157/(11.351092*min*min - 5.740650*min + 1.506631);

    // 计算βj = min{β0, β1, β2}
    real beta_j = std::min({beta0, beta1, beta2});

    // 计算δk (k=0,1,2)
    real C = 1.0;
    real delta0 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta0) ? 0 : 1;
    real delta1 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta1) ? 0 : 1;
    real delta2 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta2) ? 0 : 1;

    // 根据δ值选择插值系数
    real a_im2, a_im1, a_i, a_ip1, a_ip2;

    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a_im2 = 3.0/128.0; a_im1 = -5.0/32.0; a_i = 45.0/64.0; a_ip1 = 15.0/32.0; a_ip2 = -5.0/128.0;
    } else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a_im2 = 0.0; a_im1 = -1.0/16.0; a_i = 9.0/16.0; a_ip1 = 9.0/16.0; a_ip2 = -1.0/16.0;
    } else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
               (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a_im2 = 0.0; a_im1 = 0.0; a_i = 3.0/8.0; a_ip1 = 3.0/4.0; a_ip2 = -1.0/8.0;
    } else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a_im2 = 1.0/16.0; a_im1 = -5.0/16.0; a_i = 15.0/16.0; a_ip1 = 5.0/16.0; a_ip2 = 0.0;
    } else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a_im2 = 0.0; a_im1 = -1.0/8.0; a_i = 3.0/4.0; a_ip1 = 3.0/8.0; a_ip2 = 0.0;
    } else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a_im2 = 3.0/8.0; a_im1 = -5.0/4.0; a_i = 15.0/8.0; a_ip1 = 0.0; a_ip2 = 0.0;
    } else { // delta0==0 && delta1==0 && delta2==0
        a_im2 = 0.0; a_im1 = 0.0; a_i = 1.0; a_ip1 = 0.0; a_ip2 = 0.0;
    }

    // 计算最终插值结果
    return a_im2 * q_im2 + a_im1 * q_im1 + a_i * q_i + a_ip1 * q_ip1 + a_ip2 * q_ip2;
}



inline real whf_ai_TcnsN_myASF002_ai1(std::array<real, 5> q) {
    // 提取输入值
    real q_im2 = q[0]; // q_{i-2}
    real q_im1 = q[1]; // q_{i-1}
    real q_i   = q[2]; // q_i
    real q_ip1 = q[3]; // q_{i+1}
    real q_ip2 = q[4]; // q_{i+2}
    
    // 常数定义
    const real Cr = 0.24;
    const real xi = 1e-3;
    const real epsilon_A = (0.9 * Cr) / (1.0 - 0.9 * Cr) * xi * xi;
    const real C = 1.0;
    
    // 1. 计算局部光滑因子 β0, β1, β2
    real beta0 = std::pow(q_im2 - 2.0 * q_im1 + q_i, 2) 
               + 0.25 * std::pow(q_im2 - 4.0 * q_im1 + 3.0 * q_i, 2);
    
    real beta1 = std::pow(q_im1 - 2.0 * q_i + q_ip1, 2) 
               + 0.25 * std::pow(q_im1 - q_ip1, 2);
    
    real beta2 = std::pow(q_i - 2.0 * q_ip1 + q_ip2, 2) 
               + 0.25 * std::pow(3.0 * q_i - 4.0 * q_ip1 + q_ip2, 2);
    
    // 2. 计算全局光滑因子 τ
    real tau = std::abs(beta2 - beta0);
    
    // 3. 计算Δq值
    real delta_q_im3_2 = q_im1 - q_im2; // Δq_{i-3/2}
    real delta_q_im1_2 = q_i - q_im1;   // Δq_{i-1/2}
    real delta_q_ip1_2 = q_ip1 - q_i;   // Δq_{i+1/2}
    real delta_q_ip3_2 = q_ip2 - q_ip1; // Δq_{i+3/2}
    
    // 4. 计算η值
    real eta_im1 = (std::abs(2.0 * delta_q_im1_2 * delta_q_im3_2) + epsilon_A)
                 / (delta_q_im1_2 * delta_q_im1_2 + delta_q_im3_2 * delta_q_im3_2 + epsilon_A);
    
    real eta_i = (std::abs(2.0 * delta_q_ip1_2 * delta_q_im1_2) + epsilon_A)
               / (delta_q_ip1_2 * delta_q_ip1_2 + delta_q_im1_2 * delta_q_im1_2 + epsilon_A);
    
    real eta_ip1 = (std::abs(2.0 * delta_q_ip3_2 * delta_q_ip1_2) + epsilon_A)
                 / (delta_q_ip3_2 * delta_q_ip3_2 + delta_q_ip1_2 * delta_q_ip1_2 + epsilon_A);
    
    // 5. 计算x和辅助截断函数C_T'
    real x = std::min({1.0, eta_im1 / Cr, eta_i / Cr, eta_ip1 / Cr});
    
    real C_T_prime = 0.157 / (11.351092 * x * x - 5.740650 * x + 1.506631);
    
    // 6. 计算β_j = min(β0, β1, β2)
    real beta_j = std::min({beta0, beta1, beta2});
    
    // 7. 计算δ0, δ1, δ2
    int delta0, delta1, delta2;
    
    // 计算截断条件
    real threshold = C_T_prime * tau - C * (1.0 - C_T_prime) * beta_j;
    
    delta0 = (beta_j * tau < threshold * beta0) ? 0 : 1;
    delta1 = (beta_j * tau < threshold * beta1) ? 0 : 1;
    delta2 = (beta_j * tau < threshold * beta2) ? 0 : 1;
    
    // 8. 根据δ值选择插值系数
    real a_im2, a_im1, a_i, a_ip1, a_ip2;
    
    // 根据真值表选择系数
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        // 1,1,1
        a_im2 = 3.0/128.0;
        a_im1 = -5.0/32.0;
        a_i   = 45.0/64.0;
        a_ip1 = 15.0/32.0;
        a_ip2 = -5.0/128.0;
    }
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        // 0,1,1
        a_im2 = 0.0;
        a_im1 = -1.0/16.0;
        a_i   = 9.0/16.0;
        a_ip1 = 9.0/16.0;
        a_ip2 = -1.0/16.0;
    }
    else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
             (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        // 1,0,1 或 0,0,1
        a_im2 = 0.0;
        a_im1 = 0.0;
        a_i   = 3.0/8.0;
        a_ip1 = 3.0/4.0;
        a_ip2 = -1.0/8.0;
    }
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        // 1,1,0
        a_im2 = 1.0/16.0;
        a_im1 = -5.0/16.0;
        a_i   = 15.0/16.0;
        a_ip1 = 5.0/16.0;
        a_ip2 = 0.0;
    }
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        // 0,1,0
        a_im2 = 0.0;
        a_im1 = -1.0/8.0;
        a_i   = 3.0/4.0;
        a_ip1 = 3.0/8.0;
        a_ip2 = 0.0;
    }
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        // 1,0,0
        a_im2 = 3.0/8.0;
        a_im1 = -5.0/4.0;
        a_i   = 15.0/8.0;
        a_ip1 = 0.0;
        a_ip2 = 0.0;
    }
    else { // delta0 == 0 && delta1 == 0 && delta2 == 0
        // 0,0,0
        a_im2 = 0.0;
        a_im1 = 0.0;
        a_i   = 1.0;
        a_ip1 = 0.0;
        a_ip2 = 0.0;
    }
    
    // 9. 计算最终的插值结果
    real f_hat = a_im2 * q_im2 + a_im1 * q_im1 + a_i * q_i + a_ip1 * q_ip1 + a_ip2 * q_ip2;
    
    return f_hat;
}




//【王鸿飞】end插值格式开发（3.0）AI-新格式

// 【王鸿飞】begin插值格式终章（4.0）手搓及优化




inline real whf_TcnsN_A(std::array<real, 5> q)
{
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0/1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
              1.0/4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0/1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 
              1.0/4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0/1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 
              1.0/4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);

    // 计算自适应阈值CT_A
    real alpha1 = 10.0, alpha2 = 5.0;  //参数
    // real xi = 1e-3, C_r = 0.24;
    // real epsilon_A = (0.9*C_r)/(1.0 - 0.9*C_r)*xi*xi;
    real epsilon_A = 2.7551*1e-7;

    std::array<real, 4> delta_q;
    delta_q[0] = q[0] - q[1];
    delta_q[1] = q[1] - q[2];
    delta_q[2] = q[2] - q[3];
    delta_q[3] = q[3] - q[4];

    std::array<real, 3> eta;
    eta[0] = (std::abs(2.0 * delta_q[0] * delta_q[2]) + epsilon_A)
           / (delta_q[0] * delta_q[0] + delta_q[2] * delta_q[2] + epsilon_A);
    eta[1] = (std::abs(2.0 * delta_q[1] * delta_q[3]) + epsilon_A)
           / (delta_q[1] * delta_q[1] + delta_q[3] * delta_q[3] + epsilon_A);
    eta[2] = (std::abs(2.0 * delta_q[2] * delta_q[1]) + epsilon_A)
           / (delta_q[2] * delta_q[2] + delta_q[1] * delta_q[1] + epsilon_A);
    
    // 计算η_min值
    // real eta_min = std::min(eta[0],eta[1],eta[2]);
    real eta_min = std::min({eta[0],eta[1],eta[2]});

    // real m = 1.0 - std::min(1.0, eta_min/C_r);
    real m = 1.0 - std::min(1.0, 4.1666667 * eta_min);
    
    real g_m = (1 - m)*(1 - m)*(1 - m)*(1 - m)*(1 + 4*m);

    int beta_A = std::floor(alpha1 - alpha2*(1 - g_m));

    real CT_A;

    // 求光滑度量gamma
    const real C = 1.0;
    const real eps = 1e-40;

    std::array<real, 3> gamma;

    gamma[0] = C + tau/(beta[0] + eps);
    gamma[1] = C + tau/(beta[1] + eps);
    gamma[2] = C + tau/(beta[2] + eps);

    gamma[0] = gamma[0]*gamma[0];
    gamma[1] = gamma[1]*gamma[1];
    gamma[2] = gamma[2]*gamma[2];

    gamma[0] = gamma[0]*gamma[0]*gamma[0];
    gamma[1] = gamma[1]*gamma[1]*gamma[1];
    gamma[2] = gamma[2]*gamma[2]*gamma[2];

    real gamma_sum = gamma[0] + gamma[1] + gamma[2];
    
    switch(beta_A) {
    case 5: CT_A = 1e-5; 
    break;
    case 6: CT_A = 1e-6;
    break;
    case 7: CT_A = 1e-7;
    break;
    case 8: CT_A = 1e-8;
    break;
    case 9: CT_A = 1e-9;
    break;
    case 10: CT_A = 1e-10;
    break;
    }
    
    real rr;
    rr = CT_A*gamma_sum;
    // 构造截止函数
    unsigned short flag = 0;
    if (gamma[0] < rr) flag += 1;
    if (gamma[1] < rr) flag += 2;
    if (gamma[2] < rr) flag += 4;

    // 插值
    switch (flag) {
    case 0:  // 111
        return 3.0/128.0 * q[0] - 5.0/32.0 * q[1] + 45.0/64.0 * q[2] + 15.0/32.0 * q[3] - 5.0/128.0 * q[4];
    case 1:  // 011
        return -1.0/16.0 * q[1] + 9.0/16.0 * q[2] + 9.0/16.0 * q[3] - 1.0/16.0 * q[4];
    case 2:  // 101
    case 3:  // 001
        return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    case 4:  // 110
        return 1.0/16.0 * q[0] - 5.0/16.0 * q[1] + 15.0/16.0 * q[2] + 5.0/16.0 * q[3];
    case 5:  // 010
        return -1.0/8.0 * q[1] + 3.0/4.0 * q[2] + 3.0/8.0 * q[3];
    case 6:  // 100
        return 3.0/8.0 * q[0] - 5.0/4.0 * q[1] + 15.0/8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}


inline real whf_zyc_TcnsN_myASF002_1(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  //CT adapt begin
  // real xi = 1e-3;
  // real Cr = 0.24;
  // 此时 epsilon_A = 2.7551*1e-7 ， 1/Cr = 4.1666667

  real epsilon_A = 2.7551*1e-7;

  // 解包输入值
    real q_im2 = q[0]; // q_{i-2}
    real q_im1 = q[1]; // q_{i-1}
    real q_i   = q[2]; // q_i
    real q_ip1 = q[3]; // q_{i+1}
    real q_ip2 = q[4]; // q_{i+2}

    // 计算差分值
    real dq_im3_2 = q_im1 - q_im2; // Δq_{i-3/2}
    real dq_im1_2 = q_i - q_im1;   // Δq_{i-1/2}
    real dq_ip1_2 = q_ip1 - q_i;   // Δq_{i+1/2}
    real dq_ip3_2 = q_ip2 - q_ip1; // Δq_{i+3/2}
  // 计算η值
  real eta_im1 = (std::abs(2.0*dq_im1_2*dq_im3_2) + epsilon_A) / 
                  (std::pow(dq_im1_2, 2) + std::pow(dq_im3_2, 2) + epsilon_A);
  
  real eta_i = (std::abs(2.0*dq_ip1_2*dq_im1_2) + epsilon_A) / 
                (std::pow(dq_ip1_2, 2) + std::pow(dq_im1_2, 2) + epsilon_A);
  
  real eta_ip1 = (std::abs(2.0*dq_ip3_2*dq_ip1_2) + epsilon_A) / 
                  (std::pow(dq_ip3_2, 2) + std::pow(dq_ip1_2, 2) + epsilon_A);
  
  real eta_min = std::min({eta_im1, eta_i, eta_ip1});

  // 计算min
  // real min = std::min(1.0, eta_min / Cr);
  real min = std::min(1.0, 4.1666667 * eta_min);
  
  // 计算C_T_prime
  real C_T = 0.157/(11.351092*min*min - 5.740650*min + 1.506631);
  //CT adapt end
  
  real CT_1 = 1 - C_T;
  real tau = std::abs(beta[2] -
                      beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
  real rr = C_T * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}


inline real whf_TcnsN_COMPARE(std::array<real, 5> q)
{

    // 解包输入值ai
    const real qf_im2 = q[0]; // qf_{i-2}
    const real qf_im1 = q[1]; // qf_{i-1}
    const real qf_i   = q[2]; // qf_i
    const real qf_ip1 = q[3]; // qf_{i+1}
    const real qf_ip2 = q[4]; // qf_{i+2}

    // 局部光滑因子β my
    std::array<real, 3> beta_my;
    beta_my[0] = 1.0/1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
              1.0/4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta_my[1] = 1.0/1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 
              1.0/4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta_my[2] = 1.0/1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 
              1.0/4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 局部光滑因子β ai
    const real beta0 = std::pow(qf_im2 - 2*qf_im1 + qf_i, 2) + 
                      0.25 * std::pow(qf_im2 - 4*qf_im1 + 3*qf_i, 2);
    
    const real beta1 = std::pow(qf_im1 - 2*qf_i + qf_ip1, 2) + 
                      0.25 * std::pow(qf_im1 - qf_ip1, 2);
    
    const real beta2 = std::pow(qf_i - 2*qf_ip1 + qf_ip2, 2) + 
                      0.25 * std::pow(3*qf_i - 4*qf_ip1 + qf_ip2, 2);

    // 全局光滑因子τ my
    real tau_my = std::abs(beta_my[2] - beta_my[0]);

    // 全局光滑因子τ ai
    const real tau_ai = std::abs(beta2 - beta0);
    
    // 计算自适应阈值CT_A
    real alpha1 = 10.0, alpha2 = 5.0;  //参数
    real xi = 1e-3, C_r = 0.24;
    real epsilon_A = (0.9*C_r)/(1.0 - 0.9*C_r)*xi*xi;

    // 计算差分my
    std::array<real, 4> delta_q;
    delta_q[0] = q[0] - q[1];
    delta_q[1] = q[1] - q[2];
    delta_q[2] = q[2] - q[3];
    delta_q[3] = q[3] - q[4];

    // 计算差分ai
    const real dqf_im3_2 = qf_im1 - qf_im2; // Δqf_{i-3/2}
    const real dqf_im1_2 = qf_i - qf_im1;   // Δqf_{i-1/2}
    const real dqf_ip1_2 = qf_ip1 - qf_i;   // Δqf_{i+1/2}
    const real dqf_ip3_2 = qf_ip2 - qf_ip1; // Δqf_{i+3/2}

    // 计算η值my
    std::array<real, 3> eta_my;
    eta_my[0] = (std::abs(2.0 * delta_q[0] * delta_q[2]) + epsilon_A)
           / (delta_q[0] * delta_q[0] + delta_q[2] * delta_q[2] + epsilon_A);
    eta_my[1] = (std::abs(2.0 * delta_q[1] * delta_q[3]) + epsilon_A)
           / (delta_q[1] * delta_q[1] + delta_q[3] * delta_q[3] + epsilon_A);
    eta_my[2] = (std::abs(2.0 * delta_q[2] * delta_q[1]) + epsilon_A)
           / (delta_q[2] * delta_q[2] + delta_q[1] * delta_q[1] + epsilon_A);
    
    // 计算η值ai
    const real eta_im1 = (std::abs(2*dqf_im1_2*dqf_im3_2) + epsilon_A) /
                         (std::pow(dqf_im1_2, 2) + std::pow(dqf_im3_2, 2) + epsilon_A);
    
    const real eta_i = (std::abs(2*dqf_ip1_2*dqf_im1_2) + epsilon_A) /
                      (std::pow(dqf_ip1_2, 2) + std::pow(dqf_im1_2, 2) + epsilon_A);
    
    const real eta_ip1 = (std::abs(2*dqf_ip3_2*dqf_ip1_2) + epsilon_A) /
                         (std::pow(dqf_ip3_2, 2) + std::pow(dqf_ip1_2, 2) + epsilon_A);
           
    // my
      // 计算η_min值
    real eta_min_my = std::min({eta_my[0],eta_my[1],eta_my[2]});
      // 计算m
    real m_my = 1.0 - std::min(1.0, eta_min_my/C_r);
      // 计算gm
    real gm_my = (1 - m_my)*(1 - m_my)*(1 - m_my)*(1 - m_my)*(1 + 4*m_my);
      // 计算β_A
    int beta_A_my = std::floor(alpha1 - alpha2*(1 - gm_my));
      // 计算CT_A
    real CT_A_my;

    const real C = 1.0;
    const real eps = 1e-40;

    std::array<real, 3> gamma_my;

    gamma_my[0] = C + tau_my/(beta_my[0] + eps);
    gamma_my[1] = C + tau_my/(beta_my[1] + eps);
    gamma_my[2] = C + tau_my/(beta_my[2] + eps);

    gamma_my[0] = gamma_my[0]*gamma_my[0];
    gamma_my[1] = gamma_my[1]*gamma_my[1];
    gamma_my[2] = gamma_my[2]*gamma_my[2];

    gamma_my[0] = gamma_my[0]*gamma_my[0]*gamma_my[0];
    gamma_my[1] = gamma_my[1]*gamma_my[1]*gamma_my[1];
    gamma_my[2] = gamma_my[2]*gamma_my[2]*gamma_my[2];

    real gamma_sum_my = gamma_my[0] + gamma_my[1] + gamma_my[2];
    
    switch(beta_A_my) {
    case 5: CT_A_my = 1e-5; 
    break;
    case 6: CT_A_my = 1e-6;
    break;
    case 7: CT_A_my = 1e-7;
    break;
    case 8: CT_A_my = 1e-8;
    break;
    case 9: CT_A_my = 1e-9;
    break;
    case 10: CT_A_my = 1e-10;
    break;
    }
    
    // ai
      // 计算η_min值
    const real eta_min_ai = std::min({eta_im1, eta_i, eta_ip1});
      // 计算m
    const real m_ai = 1.0 - std::min(1.0, eta_min_ai/C_r);
      // 计算g(m)
    const real g_m_ai = std::pow(1.0 - m_ai, 4) * (1.0 + 4.0*m_ai);
      // 计算β_A
    const real beta_A_ai = alpha1 - alpha2 * (1.0 - g_m_ai);
      // 计算CT_A
    const real C_TA = std::pow(10.0, -std::floor(beta_A_ai));

      // 计算γ_k ai
    const real gamma0 = std::pow(C + tau_ai/(beta0 + eps), 6);
    const real gamma1 = std::pow(C + tau_ai/(beta1 + eps), 6);
    const real gamma2 = std::pow(C + tau_ai/(beta2 + eps), 6);

      // 归一化光滑度量χ_k ai
    const real gamma_sum_ai = gamma0 + gamma1 + gamma2;
    const real chi0 = gamma0 / gamma_sum_ai;
    const real chi1 = gamma1 / gamma_sum_ai;
    const real chi2 = gamma2 / gamma_sum_ai;



    // 求权重my
    real rr;
    rr = CT_A_my*gamma_sum_my;
    real interpolation_my = 0.0;
      // 构造截止函数my
    unsigned short flag = 0;
    if (gamma_my[0] < rr) flag += 1;
    if (gamma_my[1] < rr) flag += 2;
    if (gamma_my[2] < rr) flag += 4;

    switch (flag) {
    case 0:  // 111
        interpolation_my = 3.0/128.0 * q[0] - 5.0/32.0 * q[1] + 45.0/64.0 * q[2] + 15.0/32.0 * q[3] - 5.0/128.0 * q[4];
    case 1:  // 011
        interpolation_my = -1.0/16.0 * q[1] + 9.0/16.0 * q[2] + 9.0/16.0 * q[3] - 1.0/16.0 * q[4];
    case 2:  // 101
    case 3:  // 001
        interpolation_my = 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    case 4:  // 110
        interpolation_my = 1.0/16.0 * q[0] - 5.0/16.0 * q[1] + 15.0/16.0 * q[2] + 5.0/16.0 * q[3];
    case 5:  // 010
        interpolation_my = -1.0/8.0 * q[1] + 3.0/4.0 * q[2] + 3.0/8.0 * q[3];
    case 6:  // 100
        interpolation_my = 3.0/8.0 * q[0] - 5.0/4.0 * q[1] + 15.0/8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        interpolation_my = q[2];
    }
    

    // 求权重ai
    const int delta0 = (chi0 < C_TA) ? 0 : 1;
    const int delta1 = (chi1 < C_TA) ? 0 : 1;
    const int delta2 = (chi2 < C_TA) ? 0 : 1;

    real interpolation_ai = 0.0;

      // 根据δ值选择系数
    real a_im2, a_im1, a_i, a_ip1, a_ip2;
    
    if (delta0 == 1 && delta1 == 1 && delta2 == 1) {
        a_im2 = 3.0/128.0;
        a_im1 = -5.0/32.0;
        a_i = 45.0/64.0;
        a_ip1 = 15.0/32.0;
        a_ip2 = -5.0/128.0;
    }
    else if (delta0 == 0 && delta1 == 1 && delta2 == 1) {
        a_im2 = 0.0;
        a_im1 = -1.0/16.0;
        a_i = 9.0/16.0;
        a_ip1 = 9.0/16.0;
        a_ip2 = -1.0/16.0;
    }
    else if ((delta0 == 1 && delta1 == 0 && delta2 == 1) || 
             (delta0 == 0 && delta1 == 0 && delta2 == 1)) {
        a_im2 = 0.0;
        a_im1 = 0.0;
        a_i = 3.0/8.0;
        a_ip1 = 3.0/4.0;
        a_ip2 = -1.0/8.0;
    }
    else if (delta0 == 1 && delta1 == 1 && delta2 == 0) {
        a_im2 = 1.0/16.0;
        a_im1 = -5.0/16.0;
        a_i = 15.0/16.0;
        a_ip1 = 5.0/16.0;
        a_ip2 = 0.0;
    }
    else if (delta0 == 0 && delta1 == 1 && delta2 == 0) {
        a_im2 = 0.0;
        a_im1 = -1.0/8.0;
        a_i = 3.0/4.0;
        a_ip1 = 3.0/8.0;
        a_ip2 = 0.0;
    }
    else if (delta0 == 1 && delta1 == 0 && delta2 == 0) {
        a_im2 = 3.0/8.0;
        a_im1 = -5.0/4.0;
        a_i = 15.0/8.0;
        a_ip1 = 0.0;
        a_ip2 = 0.0;
    }
    else {  // delta0 == 0 && delta1 == 0 && delta2 == 0
        a_im2 = 0.0;
        a_im1 = 0.0;
        a_i = 1.0;
        a_ip1 = 0.0;
        a_ip2 = 0.0;
    }

    // inter
    interpolation_ai = a_im2 * qf_im2 + a_im1 * qf_im1 + a_i * qf_i + a_ip1 * qf_ip1 + a_ip2 * qf_ip2;

    // for marching
    switch (flag) {
    case 0:  // 111
        return 3.0/128.0 * q[0] - 5.0/32.0 * q[1] + 45.0/64.0 * q[2] + 15.0/32.0 * q[3] - 5.0/128.0 * q[4];
    case 1:  // 011
        return -1.0/16.0 * q[1] + 9.0/16.0 * q[2] + 9.0/16.0 * q[3] - 1.0/16.0 * q[4];
    case 2:  // 101
    case 3:  // 001
        return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    case 4:  // 110
        return 1.0/16.0 * q[0] - 5.0/16.0 * q[1] + 15.0/16.0 * q[2] + 5.0/16.0 * q[3];
    case 5:  // 010
        return -1.0/8.0 * q[1] + 3.0/4.0 * q[2] + 3.0/8.0 * q[3];
    case 6:  // 100
        return 3.0/8.0 * q[0] - 5.0/4.0 * q[1] + 15.0/8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }

}

// 【王鸿飞】end插值格式终章（4.0）手搓及优化
