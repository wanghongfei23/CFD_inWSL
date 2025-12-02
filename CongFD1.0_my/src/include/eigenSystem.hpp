/**
 * @file eigenSystem.hpp
 * @brief 定义特征系统相关的类，用于欧拉方程的特征分解
 */

#pragma once
#include<array>
#include<macro.hpp>

enum{X,Y};

/**
 * @brief 二维欧拉方程的特征系统类
 * 
 * 该类实现了二维欧拉方程的特征分解，包括将原始变量转换为特征变量，
 * 以及将特征变量转换回原始变量的功能。
 */
class eigensystemEuler2D
{
    public:
    /**
     * @brief 默认构造函数
     */
    eigensystemEuler2D(){};
    
    /**
     * @brief 构造函数，使用原始变量和法向量初始化
     * @param prim 原始变量数组 [rho, u, v, p]
     * @param norm_ 法向量数组
     */
    eigensystemEuler2D(const std::array<real,4> & prim,const std::array<real,3> & norm_);
    
    /**
     * @brief 构造函数，使用左右两侧原始变量和法向量初始化
     * @param priml 左侧原始变量数组 [rho, u, v, p]
     * @param primr 右侧原始变量数组 [rho, u, v, p]
     * @param norm_ 法向量数组
     */
    eigensystemEuler2D(const std::array<real,4> &priml,const std::array<real,4> &primr,const std::array<real,3> & norm_);
    
    /**
     * @brief 将原始变量转换为特征变量
     * @param prim 原始变量数组 [rho, u, v, p]
     * @return 特征变量数组
     */
    std::array<real,4> primToChar(const std::array<real,4> & prim);
    
    /**
     * @brief 将特征变量转换为原始变量
     * @param chars 特征变量数组
     * @return 原始变量数组 [rho, u, v, p]
     */
    std::array<real,4> charToPrim(const std::array<real,4> & chars);
    

    private:
    real r;                     ///< 密度
    real u;                     ///< x方向速度
    real v;                     ///< y方向速度
    real p;                     ///< 压力
    real gamma;                 ///< 比热比
    real ek;                    ///< 动能
    real h;                     ///< 焓
    real c;                     ///< 声速
    real Vn;                    ///< 法向速度
    bool xOrY;                  ///< 方向标志
    std::array<real,3> norm;    ///< 法向量
    std::array<real,4*4> leftEig;   ///< 左特征矩阵
    std::array<real,4*4> rightEig;  ///< 右特征矩阵
};

/**
 * @brief 一维欧拉方程的特征系统类
 * 
 * 该类实现了一维欧拉方程的特征分解，包括将原始变量转换为特征变量，
 * 以及将特征变量转换回原始变量的功能。
 */
class eigensystemEuler1D
{
    public:
    /**
     * @brief 构造函数，使用原始变量初始化
     * @param prim 原始变量数组 [rho, u, p]
     */
    eigensystemEuler1D(const std::array<real,3> & prim);
    
    /**
     * @brief 构造函数，使用左右两侧原始变量初始化
     * @param priml 左侧原始变量数组 [rho, u, p]
     * @param primr 右侧原始变量数组 [rho, u, p]
     */
    eigensystemEuler1D(const std::array<real,3> &priml,const std::array<real,3> &primr);
    
    /**
     * @brief 将原始变量转换为特征变量
     * @param prim 原始变量数组 [rho, u, p]
     * @return 特征变量数组
     */
    std::array<real,3> primToChar(const std::array<real,3> & prim);
    
    /**
     * @brief 将特征变量转换为原始变量
     * @param chars 特征变量数组
     * @return 原始变量数组 [rho, u, p]
     */
    std::array<real,3> charToPrim(const std::array<real,3> & chars);
    

    private:
    real r;                     ///< 密度
    real u;                     ///< 速度
    real p;                     ///< 压力
    real gamma;                 ///< 比热比
    real ek;                    ///< 动能
    real h;                     ///< 焓
    real c;                     ///< 声速
    bool xOrY=false;                  ///< 方向标志
    std::array<real,3*3> leftEig;   ///< 左特征矩阵
    std::array<real,3*3> rightEig;  ///< 右特征矩阵
};

// std::array<real,4> characteristicDecomposition(std::array<real,4> prim)
// {
//     return{0,0,0,0};
// }
// std::array<real,4*4> leftEigen(const std::array<real,4>& prim,const std::array<real,3>& norm)
// {
//     std::array<real,4*4> res;
//     
//     real r=prim[0],u=prim[1],v=prim[2],p=prim[3];
//     real gamma=GAMMA,ek=(u*u+v*v)/2;
//     real h=p/r*gamma/(1-gamma);
//     real Vn=norm[0]*u+norm[1]*v;
//     real c=std::sqrt(gamma*p/r);
//     //first line
//     res[0]=-norm[Y]*u+norm[X]*v;
//     res[1]=norm[Y];
//     res[2]=-norm[X];
//     res[3]=0;
//
//     //second line
//     res[4]=h-ek;
//     res[5]=u;
//     res[6]=v;
//     res[7]=-1;
//
//     //third line
//     res[8 ]=(Vn/c+ek/h)/2;
//     res[9 ]=(-norm[X]/c-u/h)/2;
//     res[10]=(-norm[Y]/c-v/h)/2;
//     res[11]=1.0/(2*h);
//
//     //fourth line
//     res[12]=(-Vn/c+ek/h)/2;
//     res[13]=(norm[X]/c-u/h)/2;
//     res[14]=(norm[Y]/c-v/h)/2;
//     res[15]=1.0/(2*h);
//
//     return res;
// }
//
// std::array<real,4*4> rightEigen(const std::array<real,4>& prim,const std::array<real,3>& norm)
// {
//     std::array<real,4*4> res;
//     enum{X,Y};
//     real r=prim[0],u=prim[1],v=prim[2],p=prim[3];
//     real gamma=GAMMA,ek=(u*u+v*v)/2;
//     real h=p/r*gamma/(1-gamma);
//     real Vn=norm[0]*u+norm[1]*v;
//     real c=std::sqrt(gamma*p/r);
//     //first line
//     res[0]=0;
//     res[1]=1/h;
//     res[2]=1;
//     res[3]=1;
//
//     //second line
//     res[4]=norm[Y];
//     res[5]=u/h;
//     res[6]=u-norm[X]*c;
//     res[7]=u+norm[X]*c;
//
//     //third line
//     res[8 ]=-norm[X];
//     res[9 ]=v/h;
//     res[10]=v-norm[Y]*c;
//     res[11]=v+norm[Y]*c;
//
//     //fourth line
//     res[12]=norm[Y]*u-norm[X]*v;
//     res[13]=ek/h;
//     res[14]=h+ek-Vn*c;
//     res[15]=h+ek+Vn*c;
//     return res;
// }