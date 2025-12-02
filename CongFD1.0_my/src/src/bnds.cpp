/**
 * @file bnds.cpp
 * @brief 边界条件管理类的实现文件
 */

#include"Bnds.hpp"

/**
 * @brief 获取指定方向的一维边界
 * @param idim 方向索引
 * @param i 第一索引
 * @param j 第二索引
 * @return 包含左右边界的数组
 */
std::array<std::shared_ptr<OneDBnd>,2> Bnds::getOneDBnd(int idim,int i,int j)
{
    std::array<std::shared_ptr<OneDBnd>,2> res;
    int index;
    if (idim==1)
    {
        index=(i+j*iMax[2])*2;
        res[0]=oneDBnds.at(index);
        res[1]=oneDBnds.at(index+1);
    }
    else if (idim==2)
    {
        index=(iMax[1]*iMax[2]+i+j*iMax[2])*2;
        res[0]=oneDBnds.at(index);
        res[1]=oneDBnds.at(index+1);
    }
    else if (idim==3)
    {
        index=(iMax[1]*iMax[2]+iMax[0]*iMax[2]+i+j*iMax[1])*2;
        res[0]=oneDBnds.at(index);
        res[1]=oneDBnds.at(index+1);
    }
    return res;
    
    
    
}

/**
 * @brief 更新所有边界条件
 */
void Bnds::update()
{
    for (auto ibnd:oneDBnds)
    {
        ibnd->update();
    }
    
}