#include "inter_scheme.hpp"
#include <algorithm>
#include <cmath>

real whf_TcnsN_A(std::array<real, 5> q)
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
    real xi = 1e-3, C_r = 0.24;
    real epsilon_A = (0.9*C_r)/(1.0 - 0.9*C_r)*xi*xi;

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

    real m = 1.0 - std::min(1.0, eta_min/C_r);
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

real whf_ai_TcnsN_A_1(std::array<real, 5> q) {
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

    // 计算最终插值
    return a_im2 * f_im2 + a_im1 * f_im1 + a_i * f_i + a_ip1 * f_ip1 + a_ip2 * f_ip2;
}


real whf_TcnsN_COMPARE(std::array<real, 5> q)
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
