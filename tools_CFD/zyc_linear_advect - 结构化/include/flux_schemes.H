#ifndef FLUX_SCHEMES_H
#define FLUX_SCHEMES_H

#include <algorithm>
#include <array>
#include <cmath>
#include <quadmath.h>

// 定义高精度浮点类型real_t为__float128
typedef __float128 real_t;

/**
 * @brief TENO5格式的实现，由CongZ提出
 * @param q 包含5个点的数组，用于重构
 * @return 重构后的界面值
 * 
 * 该函数实现了TENO5（Targeted ENO）格式，用于数值求解双曲守恒律方程
 * TENO格式通过检测间断并相应调整权重来提高精度
 */
inline real_t Teno5_CongZ(const std::array<real_t, 5>& q)
{
    // real_t eps = 1.0e-40Q;
    // 计算三个子模板的光滑指示子β
    std::array<real_t, 3> beta = {
        1.0Q / 1.0Q * powq(1.0Q * q[0] - 2.0Q * q[1] + 1.0Q * q[2], 2) + 1.0Q / 4.0Q * powq(1.0Q * q[0] - 4.0Q * q[1] + 3.0Q * q[2], 2),
        1.0Q / 1.0Q * powq(1.0Q * q[1] - 2.0Q * q[2] + 1.0Q * q[3], 2) + 1.0Q / 4.0Q * powq(1.0Q * q[1] + 0.0Q * q[2] - 1.0Q * q[3], 2),
        1.0Q / 1.0Q * powq(1.0Q * q[2] - 2.0Q * q[3] + 1.0Q * q[4], 2) + 1.0Q / 4.0Q * powq(3.0Q * q[2] - 4.0Q * q[3] + 1.0Q * q[4], 2)
    };
    // 找到最小的光滑指示子索引
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    // constexpr real_t CT = powq(1.5 * 1.0e-5, 1.0 / 6.0);
    // constexpr real_t CT_1 = 1.0 - CT;
    // 计算阈值参数
    const real_t CT = powq(1.5 * 1.0e-5, 1.0 / 6.0);
    const real_t CT_1 = 1.0 - CT;
    // 计算全局光滑指示子τ
    real_t tau = std::abs(beta[2] - beta[0]);
    // 计算判定参数
    real_t rr = CT * tau - CT_1 * beta[minBeta];
    real_t ll = tau * beta[minBeta];
    // 根据光滑性判断是否使用模板
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    // 根据flag值选择不同的重构系数
    switch (flag) {
    // switch (0) {
    case 0:
        /* 1,1,1 */
        return 3.0Q / 128.0Q * q[0] - 5.0Q / 32.0Q * q[1] + 45.0Q / 64.0Q * q[2] + 15.0Q / 32.0Q * q[3] - 5.0Q / 128.0Q * q[4];
    case 1:
        /* 0,1,1 */
        return -1.0Q / 16.0Q * q[1] + 9.0Q / 16.0Q * q[2] + 9.0Q / 16.0Q * q[3] - 1.0Q / 16.0Q * q[4];
    case 2:
        /* 1,0,1 */
        return 3.0Q / 8.0Q * q[2] + 3.0Q / 4.0Q * q[3] - 1.0Q / 8.0Q * q[4];
    case 3:
        /* 0,0,1 */
        return 3.0Q / 8.0Q * q[2] + 3.0Q / 4.0Q * q[3] - 1.0Q / 8.0Q * q[4];
    case 4:
        /* 1,1,0 */
        return 1.0Q / 16.0Q * q[0] - 5.0Q / 16.0Q * q[1] + 15.0Q / 16.0Q * q[2] + 5.0Q / 16.0Q * q[3];
    case 5:
        /* 0,1,0 */
        return -1.0Q / 8.0Q * q[1] + 3.0Q / 4.0Q * q[2] + 3.0Q / 8.0Q * q[3];
    case 6:
        /* 1,0,0 */
        return 3.0Q / 8.0Q * q[0] - 5.0Q / 4.0Q * q[1] + 15.0Q / 8.0Q * q[2];
    default:
        /* 0,0,0 */
        return q[2];
    }
}

/**
 * @brief TCNS格式的自适应切换函数实现（版本1）
 * @param q 包含5个点的数组，用于重构
 * @return 重构后的界面值
 * 
 * 该函数实现了TCNS（Targeted Compact Nonlinear Scheme）格式的自适应切换函数
 * 通过检测解的局部特征来自适应调整数值格式参数
 */
inline real_t whf_ai_TcnsN_AS_1(std::array<real_t, 5> q) {
    // 常量定义
    constexpr real_t xi = 1e-3;           // 小参数，用于避免除零
    constexpr real_t Cr = 0.24;           // 切换参数
    constexpr real_t alpha1 = 10.0;       // 自适应参数1
    constexpr real_t alpha2 = 5.0;        // 自适应参数2
    constexpr real_t qA = 6.0;            // 指数参数
    constexpr real_t C = 1.0;             // 常数
    constexpr real_t epsilon_A = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;  // 自适应参数ε_A

    // 解包输入值
    real_t q_im2 = q[0]; // q_{i-2}
    real_t q_im1 = q[1]; // q_{i-1}
    real_t q_i   = q[2]; // q_i
    real_t q_ip1 = q[3]; // q_{i+1}
    real_t q_ip2 = q[4]; // q_{i+2}

    // 计算差分值
    real_t dq_im3_2 = q_im1 - q_im2; // Δq_{i-3/2}
    real_t dq_im1_2 = q_i - q_im1;   // Δq_{i-1/2}
    real_t dq_ip1_2 = q_ip1 - q_i;   // Δq_{i+1/2}
    real_t dq_ip3_2 = q_ip2 - q_ip1; // Δq_{i+3/2}

    // 计算局部光滑因子
    real_t beta0 = powq(q_im2 - 2*q_im1 + q_i, 2) + 
                 0.25 * powq(q_im2 - 4*q_im1 + 3*q_i, 2);
    
    real_t beta1 = powq(q_im1 - 2*q_i + q_ip1, 2) + 
                 0.25 * powq(q_im1 - q_ip1, 2);
    
    real_t beta2 = powq(q_i - 2*q_ip1 + q_ip2, 2) + 
                 0.25 * powq(3*q_i - 4*q_ip1 + q_ip2, 2);

    // 计算全局光滑因子
    real_t tau = std::abs(beta2 - beta0);

    // 计算η值
    real_t eta_im1 = (std::abs(2*dq_im1_2*dq_im3_2) + epsilon_A) / 
                   (powq(dq_im1_2, 2) + powq(dq_im3_2, 2) + epsilon_A);
    
    real_t eta_i = (std::abs(2*dq_ip1_2*dq_im1_2) + epsilon_A) / 
                 (powq(dq_ip1_2, 2) + powq(dq_im1_2, 2) + epsilon_A);
    
    real_t eta_ip1 = (std::abs(2*dq_ip3_2*dq_ip1_2) + epsilon_A) / 
                   (powq(dq_ip3_2, 2) + powq(dq_ip1_2, 2) + epsilon_A);
    
    real_t eta_min = std::min({eta_im1, eta_i, eta_ip1});

    // 计算m和g(m)
    real_t m = 1.0Q - std::min(1.0Q, eta_min / Cr);
    real_t g_m = powq(1.0Q - m, 4) * (1.0Q + 4.0Q*m);
    
    // 计算自适应参数
    real_t beta_A = alpha1 - alpha2 * (1.0Q - g_m);
    real_t C_TA = powq(10.0, -floorq(beta_A));
    real_t C_T_prime = powq(1.5 * C_TA / (1.0 - C_TA), 1.0/qA);

    // 确定β_j
    real_t beta_j = std::min({beta0, beta1, beta2});

    // 计算δ_k
    int delta0 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta0) ? 0 : 1;
    int delta1 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta1) ? 0 : 1;
    int delta2 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta2) ? 0 : 1;

    // 根据δ值选择插值系数
    real_t a0, a1, a2, a3, a4;
    
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

/**
 * @brief TCNS格式的自适应切换函数实现（myASF002版本1）
 * @param q 包含5个点的数组，用于重构
 * @return 重构后的界面值
 * 
 * 这是TCNS格式的另一个自适应切换函数实现版本
 * 使用不同的自适应参数计算方法
 */
inline real_t whf_ai_TcnsN_myASF002_1(std::array<real_t, 5> q) {
    // 常量定义
    constexpr real_t xi = 1e-3;           // 小参数
    constexpr real_t Cr = 0.24;           // 切换参数
    constexpr real_t C = 1.0;             // 常数
    constexpr real_t epsilon_A = (0.9 * Cr) / (1 - 0.9 * Cr) * xi * xi;  // 自适应参数ε_A

    // 解包输入值
    real_t q_im2 = q[0]; // q_{i-2}
    real_t q_im1 = q[1]; // q_{i-1}
    real_t q_i   = q[2]; // q_i
    real_t q_ip1 = q[3]; // q_{i+1}
    real_t q_ip2 = q[4]; // q_{i+2}

    // 计算差分值
    real_t dq_im3_2 = q_im1 - q_im2; // Δq_{i-3/2}
    real_t dq_im1_2 = q_i - q_im1;   // Δq_{i-1/2}
    real_t dq_ip1_2 = q_ip1 - q_i;   // Δq_{i+1/2}
    real_t dq_ip3_2 = q_ip2 - q_ip1; // Δq_{i+3/2}

    // 计算局部光滑因子
    real_t beta0 = powq(q_im2 - 2*q_im1 + q_i, 2) + 
                 0.25 * powq(q_im2 - 4*q_im1 + 3*q_i, 2);
    
    real_t beta1 = powq(q_im1 - 2*q_i + q_ip1, 2) + 
                 0.25 * powq(q_im1 - q_ip1, 2);
    
    real_t beta2 = powq(q_i - 2*q_ip1 + q_ip2, 2) + 
                 0.25 * powq(3*q_i - 4*q_ip1 + q_ip2, 2);

    // 计算全局光滑因子
    real_t tau = std::abs(beta2 - beta0);

    // 计算η值
    real_t eta_im1 = (std::abs(2*dq_im1_2*dq_im3_2) + epsilon_A) / 
                   (powq(dq_im1_2, 2) + powq(dq_im3_2, 2) + epsilon_A);
    
    real_t eta_i = (std::abs(2*dq_ip1_2*dq_im1_2) + epsilon_A) / 
                 (powq(dq_ip1_2, 2) + powq(dq_im1_2, 2) + epsilon_A);
    
    real_t eta_ip1 = (std::abs(2*dq_ip3_2*dq_ip1_2) + epsilon_A) / 
                   (powq(dq_ip3_2, 2) + powq(dq_ip1_2, 2) + epsilon_A);
    
    real_t eta_min = std::min({eta_im1, eta_i, eta_ip1});

    // 计算min参数
    real_t min = std::min(1.0Q, eta_min / Cr);
    
    // 计算C_T_prime（自适应阈值参数）
    real_t C_T_prime = 0.157/(11.351092*min*min - 5.740650*min + 1.506631);

    // 确定β_j
    real_t beta_j = std::min({beta0, beta1, beta2});

    // 计算δ_k
    int delta0 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta0) ? 0 : 1;
    int delta1 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta1) ? 0 : 1;
    int delta2 = (beta_j * tau < (C_T_prime * tau - C * (1 - C_T_prime) * beta_j) * beta2) ? 0 : 1;

    // 根据δ值选择插值系数
    real_t a0, a1, a2, a3, a4;
    
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

/**
 * @brief TCNS格式的A版本实现
 * @param q 包含5个点的数组，用于重构
 * @return 重构后的界面值
 * 
 * 该函数实现了TCNS（Targeted Compact Nonlinear Scheme）格式的A版本
 * 通过自适应阈值参数来调整数值格式
 */
inline real_t whf_TcnsN_A(std::array<real_t, 5> q)
{
    // 局部光滑因子β
    std::array<real_t, 3> beta;
    beta[0] = 1.0/1.0 * powq(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
              1.0/4.0 * powq(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0/1.0 * powq(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 
              1.0/4.0 * powq(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0/1.0 * powq(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 
              1.0/4.0 * powq(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real_t tau = std::abs(beta[2] - beta[0]);

    // 计算自适应阈值CT_A
    real_t alpha1 = 10.0, alpha2 = 5.0;  //自适应参数
    real_t xi = 1e-3, C_r = 0.24;        //小参数和切换参数
    real_t epsilon_A = (0.9*C_r)/(1.0 - 0.9*C_r)*xi*xi;  //自适应参数ε_A

    // 计算差分值
    std::array<real_t, 4> delta_q;
    delta_q[0] = q[0] - q[1];
    delta_q[1] = q[1] - q[2];
    delta_q[2] = q[2] - q[3];
    delta_q[3] = q[3] - q[4];

    // 计算η值（用于检测解的局部特征）
    std::array<real_t, 3> eta;
    eta[0] = (std::abs(2.0 * delta_q[0] * delta_q[2]) + epsilon_A)
           / (delta_q[0] * delta_q[0] + delta_q[2] * delta_q[2] + epsilon_A);
    eta[1] = (std::abs(2.0 * delta_q[1] * delta_q[3]) + epsilon_A)
           / (delta_q[1] * delta_q[1] + delta_q[3] * delta_q[3] + epsilon_A);
    eta[2] = (std::abs(2.0 * delta_q[2] * delta_q[1]) + epsilon_A)
           / (delta_q[2] * delta_q[2] + delta_q[1] * delta_q[1] + epsilon_A);
    
    // 计算η_min值（最小的η值）
    real_t eta_min = std::min({eta[0],eta[1],eta[2]});

    // 计算自适应函数参数
    real_t m = 1.0 - std::min(1.0Q, eta_min/C_r);
    real_t g_m = (1 - m)*(1 - m)*(1 - m)*(1 - m)*(1 + 4*m);

    // 计算beta_A参数
    int beta_A = (int)floorq(alpha1 - alpha2*(1 - g_m));

    // 自适应阈值参数CT_A
    real_t CT_A;

    // 求光滑度量gamma（用于构造截止函数）
    const real_t C = 1.0;
    const real_t eps = 1e-40;

    std::array<real_t, 3> gamma;

    gamma[0] = C + tau/(beta[0] + eps);
    gamma[1] = C + tau/(beta[1] + eps);
    gamma[2] = C + tau/(beta[2] + eps);

    // 对gamma进行幂运算处理
    gamma[0] = powq(gamma[0], 2);
    gamma[1] = powq(gamma[1], 2);
    gamma[2] = powq(gamma[2], 2);

    gamma[0] = powq(gamma[0], 3);
    gamma[1] = powq(gamma[1], 3);
    gamma[2] = powq(gamma[2], 3);

    // 计算gamma总和
    real_t gamma_sum = gamma[0] + gamma[1] + gamma[2];
    
    // 根据beta_A值确定CT_A
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
    
    // 计算判定参数rr
    real_t rr;
    rr = CT_A*gamma_sum;
    // 构造截止函数
    unsigned short flag = 0;
    if (gamma[0] < rr) flag += 1;
    if (gamma[1] < rr) flag += 2;
    if (gamma[2] < rr) flag += 4;

    // 根据flag值进行插值
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

#endif // FLUX_SCHEMES_H