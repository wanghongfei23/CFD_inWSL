#ifndef ALGEBRAIC_FUNCTION_HPP
#define ALGEBRAIC_FUNCTION_HPP

#include <cmath>

/**
 * @brief 代数函数类，用于CFD中的插值格式测试
 * 定义函数 g_k(x) = x^2 * exp(0.75*(x-1))
 */
class AlgebraicFunction {
public:
    /**
     * @brief 计算代数函数值 g_k(x) = x^2 * exp(0.75*(x-1))
     * @param x 自变量
     * @return 函数值
     */
    static double g_k(double x);

    /**
     * @brief 计算函数的一阶导数
     * @param x 自变量
     * @return 导数值
     */
    static double g_k_derivative(double x);

    /**
     * @brief 计算函数的二阶导数
     * @param x 自变量
     * @return 二阶导数值
     */
    static double g_k_second_derivative(double x);
};

#endif // ALGEBRAIC_FUNCTION_HPP