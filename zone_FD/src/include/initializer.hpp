#pragma once
#include"block.hpp"
#include"info.hpp"
#include"equation.hpp"
#include"bnds.hpp"
#include"sp_distributor.hpp"

// 初始化器类，用于初始化计算域、方程、边界条件等
class Initializer
{
    public:
    // 默认构造函数
    Initializer();
    
    // 带Info参数的构造函数
    Initializer(Info*);
    
    // 初始化解数据
    void solInit(Block*,Data*);
    
    // 初始化均匀网格块
    void initUniformBlock(Block*);
    
    // 初始化方程求解器
    void initEqution(Equation*,Block*);

    // 初始化边界条件
    void initBnds(Bnds* bnds,Equation*,std::array<int,3> iMax,Block* block);
    
    // 初始化双马赫反射问题的边界条件
    void initDoubleMachBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block);

    // 初始化空间分布器
    void initSpDistributor(SpDistributor*,Equation*,Block*,Bnds*);

    private:
    // Info对象指针，包含计算相关信息
    Info* info;
};