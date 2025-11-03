#pragma once

#include"data.hpp"

// 方程类，用于处理守恒变量和基本变量之间的转换
class Equation
{
    public:
    // 构造函数
    Equation(){};
    
    // 变量转换函数
    void consToPrim();
    
    // 数据获取函数
    Data* getCons();
    Data* getPrim();
    Data* getRhs();

    protected:
    // 允许初始化器访问保护成员
    friend class Initializer;

    // 变量转换实现函数
    void consToPrimEuler1D();
    void consToPrimEuler2D();
    
    // 方程相关数据
    bool inited=false;        // 初始化标志
    int n,nCons,nPrim;        // 网格点数、守恒变量数、基本变量数
    int dim;                  // 问题维度
    Data* cons,*rhs,*prim;    // 守恒变量、右手端、基本变量数据指针
    EquationType type;        // 方程类型
};