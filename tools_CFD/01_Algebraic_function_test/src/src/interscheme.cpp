#include "../include/interscheme.hpp"
#include <cmath>
#include <algorithm>

// 使用命名空间简化代码
using namespace std;

/**
 * @brief WHF 自定义高阶插值格式（CT自适应权重）
 * 结合了WENO思想与自适应控制参数C_T，用于提高间断分辨率
 * 
 * @param q 输入的5个网格点值
 * @return 插值结果
 */
real whf_TCNS_AS_myF203_NoS(std::array<real, 5> q) {
    real eps = 1e-40;
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
            1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
            1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
            1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };

    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();

    // CT 自适应参数设置
    real epsilon_A = 2.7551e-7;

    std::array<real, 4> delta_q;
    delta_q[0] = q[0] - q[1];
    delta_q[1] = q[1] - q[2];
    delta_q[2] = q[2] - q[3];
    delta_q[3] = q[3] - q[4];

    real eta_im1 = (abs(2.0*delta_q[1]*delta_q[0]) + epsilon_A) / 
                   (pow(delta_q[1], 2) + pow(delta_q[0], 2) + epsilon_A);

    real eta_i = (abs(2.0*delta_q[2]*delta_q[1]) + epsilon_A) / 
                 (pow(delta_q[2], 2) + pow(delta_q[1], 2) + epsilon_A);

    real eta_ip1 = (abs(2.0*delta_q[3]*delta_q[2]) + epsilon_A) / 
                   (pow(delta_q[3], 2) + pow(delta_q[2], 2) + epsilon_A);

    real eta_min = std::min({eta_im1, eta_i, eta_ip1});
    real min_val = std::min(0.24, eta_min);
    real m = 1 - min_val/0.24;

    real mm = m*m;
    real mmm = m*mm;
    real C_T = -0.452662*mmm + 0.590811*mm - 0.012189*m + 0.022010;

    real CT_1 = 1 - C_T;
    real tau = abs(beta[2] - beta[0]);
    real rr = C_T * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];

    unsigned short flag = 0;
    if (ll < rr * beta[0]) flag += 1;
    if (ll < rr * beta[1]) flag += 2;
    if (ll < rr * beta[2]) flag += 4;

    switch (flag) {
    case 0: return 3.0/128*q[0] - 5.0/32*q[1] + 45.0/64*q[2] + 15.0/32*q[3] - 5.0/128*q[4];
    case 1: return -1.0/16*q[1] + 9.0/16*q[2] + 9.0/16*q[3] - 1.0/16*q[4];
    case 2: return 3.0/8*q[2] + 3.0/4*q[3] - 1.0/8*q[4];
    case 3: return 3.0/8*q[2] + 3.0/4*q[3] - 1.0/8*q[4];
    case 4: return 1.0/16*q[0] - 5.0/16*q[1] + 15.0/16*q[2] + 5.0/16*q[3];
    case 5: return -1.0/8*q[1] + 3.0/4*q[2] + 3.0/8*q[3];
    case 6: return 3.0/8*q[0] - 5.0/4*q[1] + 15.0/8*q[2];
    default:return q[2];
    }
}

/**
 * @brief 经典五阶WENO格式（Jiang-Shu型权重）
 * 使用非线性凸组合构造高阶精度逼近
 * 
 * @param q 输入的5个网格点值
 * @return 插值结果
 */
real weno5_JSchen(std::array<real, 5> q) {
    constexpr real eps = 1e-6;
    std::array<real, 3> gamma = {1.0/16.0, 5.0/8.0, 5.0/16.0};
    std::array<real, 3> beta;

    beta[0] = pow(q[0] - 2*q[1] + q[2], 2) + 0.25*pow(q[0] - 4*q[1] + 3*q[2], 2);
    beta[1] = pow(q[1] - 2*q[2] + q[3], 2) + 0.25*pow(q[1] - q[3], 2);
    beta[2] = pow(q[2] - 2*q[3] + q[4], 2) + 0.25*pow(3*q[2] - 4*q[3] + q[4], 2);

    std::array<real, 3> u = {
        3.0/8*q[0] - 5.0/4*q[1] + 15.0/8*q[2],
        -1.0/8*q[1] + 3.0/4*q[2] + 3.0/8*q[3],
        3.0/8*q[2] + 3.0/4*q[3] - 1.0/8*q[4]
    };

    real sumw = 0.0, result = 0.0;
    for (int i = 0; i < 3; ++i) {
        real alpha = gamma[i] / (eps + beta[i]*beta[i]);
        result += alpha * u[i];
        sumw += alpha;
    }
    return result / sumw;
}

/**
 * @brief Nicest5: 多级平滑因子最优选择五阶插值
 * 先判断最优模板长度（3/4/5点），再返回对应重构值
 * 
 * @param q 输入的5个网格点值
 * @return 插值结果
 */
real Nicest5(std::array<real, 5> q) {
    constexpr real eps = 1e-40;

    const real beta_31 = pow(q[0] - 2*q[1] + q[2], 2) + 0.25*pow(q[0] - 4*q[1] + 3*q[2], 2) + eps;
    const real beta_32 = pow(q[1] - 2*q[2] + q[3], 2) + 0.25*pow(q[1] - q[3], 2) + eps;
    const real beta_33 = pow(q[2] - 2*q[3] + q[4], 2) + 0.25*pow(3*q[2] - 4*q[3] + q[4], 2) + eps;
    const real beta_51 = 0.06 * (pow(q[0] - 4*q[1] + 6*q[2] - 4*q[3] + q[4], 2) + pow(-0.5*q[0] + q[1] - q[3] + 0.5*q[4], 2)) + eps;

    const real min_beta_s1 = std::min({beta_31, beta_32, beta_33, beta_51});

    if (min_beta_s1 == beta_51) {
        return 3.0/128*q[0] - 5.0/32*q[1] + 45.0/64*q[2] + 15.0/32*q[3] - 5.0/128*q[4];
    }

    const real beta_41 = 0.06 * (pow(q[1] - 2*q[2] + q[3], 2) + pow(-q[0] + 3*q[1] - 3*q[2] + q[3], 2)) + eps;
    const real beta_42 = 0.06 * (pow(q[2] - 2*q[3] + q[4], 2) + pow(-q[1] + 3*q[2] - 3*q[3] + q[4], 2)) + eps;

    const real min_beta_s2 = std::min({beta_41, beta_42, min_beta_s1});

    if (min_beta_s2 == beta_41) {
        return 1.0/16*q[0] - 5.0/16*q[1] + 15.0/16*q[2] + 5.0/16*q[3];
    }
    if (min_beta_s2 == beta_42) {
        return -1.0/16*q[1] + 9.0/16*q[2] + 9.0/16*q[3] - 1.0/16*q[4];
    }
    if (min_beta_s2 == beta_31) {
        return 3.0/8*q[0] - 5.0/4*q[1] + 15.0/8*q[2];
    }
    if (min_beta_s2 == beta_32) {
        return -1.0/8*q[1] + 3.0/4*q[2] + 3.0/8*q[3];
    }
    if (min_beta_s2 == beta_33) {
        return 3.0/8*q[2] + 3.0/4*q[3] - 1.0/8*q[4];
    }

    return q[2]; // fallback
}

/**
 * @brief 二阶中心差分插值格式
 * 使用最简单的二阶中心差分格式进行插值
 * 
 * @param q 输入的5个网格点值
 * @return 插值结果
 */
real secondOrderCentered(std::array<real, 5> q) {
    // 使用第二和第四个点进行二阶中心差分插值
    // 在第i+1/2处插值，使用第i和第i+1个点
    // 这里我们使用q[1]和q[3]作为左右两个点
    return 0.5 * (q[1] + q[3]);
}
