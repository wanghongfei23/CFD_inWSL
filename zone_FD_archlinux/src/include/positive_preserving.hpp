#pragma once

#include "macro.hpp"
#include <concepts>

const real EPSILON = 1e-10;

template <typename F, size_t ncomp>
concept PositiveQuantityFunc = requires(F f, const std::array<real, ncomp>& W) {
    { f(W) } -> std::same_as<real>; // 要求 f(W) 返回 real 类型
};

template <size_t ncomp, PositiveQuantityFunc<ncomp> PositiveQuantityFuncType>
 void variable_positive_limiter(
    std::array<real,ncomp>const& W, // 原始变量值
    std::array<real, ncomp>& W_interpolated, // 高阶插值后的变量值
    PositiveQuantityFuncType computePositiveQuantity, // 计算正值的函数
    real EPSILON = 1e-10 // 正值的阈值
) {
    // 计算插值后变量中的正值
    real positiveQuantity_W_interpolated = computePositiveQuantity(W_interpolated);

    // 如果插值后的正值大于或等于阈值，则无需调整
    if (positiveQuantity_W_interpolated >= EPSILON) {
        return; // 直接返回，无需调整
    }

    // 计算原始变量中的正值
    real positiveQuantity_W = computePositiveQuantity(W);

    // 如果原始变量中的正值小于阈值，则完全使用原始值
    if (positiveQuantity_W < EPSILON) {
        for (size_t i = 0; i < ncomp; ++i) {
            W_interpolated[i] = W[i];
        }
        return; // 直接返回，无需进一步计算
    }

    // 如果插值后的正值小于阈值，则计算调整系数 theta
    real theta = (EPSILON - positiveQuantity_W) / (positiveQuantity_W_interpolated - positiveQuantity_W);

    // 调整插值后的变量值
    for (size_t i = 0; i < ncomp; ++i) {
        W_interpolated[i] = (1 - theta) * W[i] + theta * W_interpolated[i];
    }
}


// 通量保正限制器函数
// template <std::size_t ncomp, PositiveQuantityFunction<ncomp> Func>
// void flux_positivity_limiter(
//     const std::array<double, ncomp>& WHLLC_plus,  // HLLC+ 状态
//     const std::array<double, ncomp>& Wstar_plus,  // **+ 状态
//     const std::array<double, ncomp>& WHLLC_minus, // HLLC- 状态
//     const std::array<double, ncomp>& Wstar_minus, // **- 状态
//     const std::array<double, ncomp>& GHLLC_plus,  // HLLC+ 通量
//     std::array<double, ncomp>& Gstar_plus,        // **+ 通量（直接修改）
//     const std::array<double, ncomp>& GHLLC_minus, // HLLC- 通量
//     std::array<double, ncomp>& Gstar_minus,       // **- 通量（直接修改）
//     Func alpha_k,                                 // 计算保正量的函数
//     double epsilon_alpha=1e-10                          // 保正量阈值
// ) {
//     // 提前计算 alpha_k 的值
//     const double alpha_HLLC_plus = alpha_k(WHLLC_plus);
//     const double alpha_star_plus = alpha_k(Wstar_plus);
//     const double alpha_HLLC_minus = alpha_k(WHLLC_minus);
//     const double alpha_star_minus = alpha_k(Wstar_minus);

//     // 检查是否需要直接跳出
//     if (alpha_HLLC_plus >= epsilon_alpha && alpha_star_plus >= epsilon_alpha &&
//         alpha_HLLC_minus >= epsilon_alpha && alpha_star_minus >= epsilon_alpha) {
//         return; // 无需修正，直接返回
//     }

//     // 计算 theta_plus
//     double theta_plus = 1.0;
//     if (alpha_HLLC_plus < epsilon_alpha) {
//         theta_plus = 0.0;
//     } else if (alpha_star_plus < epsilon_alpha) {
//         theta_plus = (epsilon_alpha - alpha_HLLC_plus) / (alpha_star_plus - alpha_HLLC_plus);
//         theta_plus = std::clamp(theta_plus, 0.0, 1.0);
//     }

//     // 计算 theta_minus
//     double theta_minus = 1.0;
//     if (alpha_HLLC_minus < epsilon_alpha) {
//         theta_minus = 0.0;
//     } else if (alpha_star_minus < epsilon_alpha) {
//         theta_minus = (epsilon_alpha - alpha_HLLC_minus) / (alpha_star_minus - alpha_HLLC_minus);
//         theta_minus = std::clamp(theta_minus, 0.0, 1.0);
//     }

//     // 取 theta 的最小值
//     double theta = std::min(theta_plus, theta_minus);

//     // 凸组合计算最终通量
//     for (std::size_t i = 0; i < ncomp; ++i) {
//         Gstar_plus[i] = (1.0 - theta) * GHLLC_plus[i] + theta * Gstar_plus[i];
//         Gstar_minus[i] = (1.0 - theta) * GHLLC_minus[i] + theta * Gstar_minus[i];
//     }
// }

