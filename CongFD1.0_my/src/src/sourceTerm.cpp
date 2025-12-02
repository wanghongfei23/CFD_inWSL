/**
 * @file sourceTerm.cpp
 * @brief 源项类的实现文件
 */

#include "SourceTerm.hpp"

/**
 * @brief 计算重力源项
 * 
 * 该函数专门用于计算Rayleigh-Taylor不稳定性问题中的重力源项。
 */
void SourceTerm::calGravitySource()
{
    //here it is still only for R-T instability case
    if(nprim!=4) 
    {
        std::cout<<"SourceTerm error: nprim incorrect\n";
    }
    real r,u,v,p;
    for (int i = 0; i < n; i++)
    {
        r=(*prim)(i,0);
        u=(*prim)(i,1);
        v=(*prim)(i,2);
        p=(*prim)(i,3);
        (*rhs)(i,2)-=r;
        (*rhs)(i,3)-=r*v;
    }
    
}

/**
 * @brief 空操作函数
 * 
 * 当没有源项需要计算时调用此函数。
 */
void SourceTerm::nothingHappened()
{

}


/**
 * @brief 构造函数
 * @param prim_ 原始变量数据指针
 * @param rhs_ 右手端数据指针
 * @param info_ Info对象指针
 */
SourceTerm::SourceTerm(Data* prim_,Data* rhs_,Info* info_)
{
    prim=prim_;
    rhs=rhs_;
    info=info_;

    auto iMax=info->icMax();
    n=iMax[0]*iMax[1]*iMax[2];
    nprim=info->nPrim();
    nCons=info->nCons();
    switch (info->sourceType)
    {
    case GRAVITY:
        calSourceMethod=(&SourceTerm::calGravitySource);
        break;
    
    default:
        calSourceMethod=(&SourceTerm::nothingHappened);
        break;
    }
}

/**
 * @brief 计算源项
 * 
 * 根据初始化时设定的源项类型，调用相应的源项计算函数。
 */
void SourceTerm::calSource()
{
    (this->*calSourceMethod)();
}