#pragma once

#include"data.hpp"

/**
 * @brief 方程类，用于处理守恒变量和基本变量之间的转换
 * 
 * Equation类提供了在守恒变量（conserved variables）和基本变量（primitive variables）
 * 之间的相互转换功能，这是计算流体力学中常见的操作。该类支持多种类型的方程，包括线性对流、
 * Burgers方程和Euler方程等。
 */
class Equation
{
    public:
    /**
     * @brief 构造函数
     */
    Equation(){};
    
    /**
     * @brief 变量转换函数
     * 
     * 根据方程类型调用相应的变量转换方法，将守恒变量转换为基本变量
     */
    void consToPrim();
    
    /**
     * @brief 获取守恒变量数据指针
     * @return 指向守恒变量数据的指针
     */
    Data* getCons();
    
    /**
     * @brief 获取基本变量数据指针
     * @return 指向基本变量数据的指针
     */
    Data* getPrim();
    
    /**
     * @brief 获取右手端(RHS)数据指针
     * @return 指向右手端数据的指针
     */
    Data* getRhs();

    protected:
    /// @cond
    // 允许初始化器访问保护成员
    friend class Initializer;
    /// @endcond

    /**
     * @brief 一维Euler方程变量转换实现函数
     * 
     * 将一维Euler方程的守恒变量转换为基本变量
     */
    void consToPrimEuler1D();
    
    /**
     * @brief 二维Euler方程变量转换实现函数
     * 
     * 将二维Euler方程的守恒变量转换为基本变量
     */
    void consToPrimEuler2D();
    
    // 方程相关数据
    bool inited=false;        ///< 初始化标志
    int n,nCons,nPrim;        ///< 网格点数、守恒变量数、基本变量数
    int dim;                  ///< 问题维度
    Data* cons,*rhs,*prim;    ///< 守恒变量、右手端、基本变量数据指针
    EquationType type;        ///< 方程类型
};