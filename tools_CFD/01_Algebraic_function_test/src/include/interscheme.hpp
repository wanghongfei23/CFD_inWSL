#ifndef INTERSCHEME_HPP
#define INTERSCHEME_HPP

#include <array>
#include <cmath>
#include <algorithm>

#define real double

// 定义插值函数指针类型
typedef real (*InterpFunc)(std::array<real, 5>);

/**
 * @brief 高精度有限差分插值格式
 * 实现了一个基于WENO思想的高阶插值格式，用于CFD计算中的数值逼近
 * 
 * @param q 包含5个点的模板值数组
 * @return 插值结果
 */
real whf_TCNS_AS_myF203_NoS(std::array<real, 5> q);

/**
 * @brief WENO5_JSchen插值格式
 * 实现了经典的五阶WENO-JS格式
 * 
 * @param q 包含5个点的模板值数组
 * @return 插值结果
 */
real weno5_JSchen(std::array<real, 5> q);

/**
 * @brief Nicest5插值格式
 * 实现了NICEST五阶格式
 * 
 * @param q 包含5个点的模板值数组
 * @return 插值结果
 */
real Nicest5(std::array<real, 5> q);

/**
 * @brief 二阶中心差分插值格式
 * 实现了简单的二阶中心差分格式
 * 
 * @param q 包含5个点的模板值数组
 * @return 插值结果
 */
real secondOrderCentered(std::array<real, 5> q);

#endif // INTERSCHEME_HPP