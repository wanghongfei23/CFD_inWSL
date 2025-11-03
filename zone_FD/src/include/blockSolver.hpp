#pragma once
#include "SourceTerm.hpp"
#include "cgnsio.hpp"
#include "initializer.hpp"

// 块求解器类，用于控制整个求解过程
class BlockSolver {
public:
    // 构造和析构函数
    BlockSolver();
    BlockSolver(Info*);
    ~BlockSolver();
    
    // 求解控制函数
    void solve(real);
    
    // 输出控制函数
    void outputGrid();
    void outputCons();
    void outputPrim();
    
    // 时间步进循环控制函数
    void stepsLoop();
    void stepsLoopCFL();
    void stepsLoopDTS();
    
    // 测试函数
    void Test();
    
    // 时间步计数器
    long timesteps = 0;

private:
    // IO处理对象
    CgnsIO cgnsIO;
    
    // 配置和数据对象
    Info* info;
    Block* block;
    Initializer* initer;
    Equation* eqn;
    Bnds* bnds;
    SpDistributor* spDis;
    Data *cons, *rhs;
    SourceTerm* sourceTerm;

    // 时间积分方法
    void RK3_SSP(real);
    void RK4_SSP(real);
    void DTS_Euler(real);
    
    // CFL计算相关函数
    std::vector<real> calLocalCFL();
    real getTimeIntervalExplicit();
    
    // 时间积分方法选择
    TimeMethod timeMethod = RK3SSP;
};