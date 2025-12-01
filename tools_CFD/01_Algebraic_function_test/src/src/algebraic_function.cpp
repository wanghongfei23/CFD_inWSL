#include "../include/algebraic_function.hpp"

/**
 * @brief 计算代数函数值 g_k(x) = x^2 * exp(0.75*(x-1))
 * @param x 自变量
 * @return 函数值
 */
double AlgebraicFunction::g_k(double x) {
    return x * x * std::exp(0.75 * (x - 1.0));
}

/**
 * @brief 计算函数的一阶导数 g_k'(x) = (2x + 0.75x^2) * exp(0.75*(x-1))
 * @param x 自变量
 * @return 导数值
 */
double AlgebraicFunction::g_k_derivative(double x) {
    return (2.0 * x + 0.75 * x * x) * std::exp(0.75 * (x - 1.0));
}

/**
 * @brief 计算函数的二阶导数
 * @param x 自变量
 * @return 二阶导数值
 */
double AlgebraicFunction::g_k_second_derivative(double x) {
    return (2.0 + 1.5 * x + 0.5625 * x * x) * std::exp(0.75 * (x - 1.0));
}