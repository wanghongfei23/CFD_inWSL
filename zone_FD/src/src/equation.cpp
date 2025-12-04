#include "equation.hpp"

/**
 * @brief 变量转换函数的实现
 * 
 * 根据方程类型调用相应的变量转换方法，将守恒变量转换为基本变量
 */
void Equation::consToPrim()
{
    switch (type)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        std::copy(cons->begin(),cons->end(),prim->begin());
        break;
    case EULER:
        {
            if (dim==1) consToPrimEuler1D();
            else if(dim==2) 
            {
                consToPrimEuler2D();
            }
        }
        break;
    
    default:
        break;
    }
}

/**
 * @brief 一维Euler方程变量转换的具体实现
 * 
 * 将一维Euler方程的守恒变量(r, ru, rE)转换为基本变量(r, u, p)
 * 其中r是密度，u是速度，p是压力，E是总能
 */
void Equation::consToPrimEuler1D()
{
    if(nCons!=3,nPrim!=3)
    {
        std::cout<<"Equation error: Euler 1d equation variable number error \n";
    }

    for (int i = 0; i < n; i++)
    {
        real r=(*cons)(i,0);
        real ru=(*cons)(i,1);
        real rE=(*cons)(i,2);
        real u=ru/r;
        real E=rE/r;
        real e=-u*u/2+E;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=p;
    }
}

/**
 * @brief 二维Euler方程变量转换的具体实现
 * 
 * 将二维Euler方程的守恒变量(r, ru, rv, rE)转换为基本变量(r, u, v, p)
 * 其中r是密度，u和v是x和y方向的速度，p是压力，E是总能
 */
void Equation::consToPrimEuler2D()
{
    if(nCons!=3,nPrim!=4)
    {
        std::cout<<"Equation error: Euler 1d equation variable number error \n";
    }

    for (int i = 0; i < n; i++)
    {
        real r=(*cons)(i,0);
        real ru=(*cons)(i,1);
        real rv=(*cons)(i,2);
        real rE=(*cons)(i,3);
        real u=ru/r;
        real v=rv/r;
        real E=rE/r;
        real q2=(u*u+v*v)/2;
        real e=-q2+E;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=v;
        (*prim)(i,3)=p;
    }
}


/**
 * @brief 获取基本变量数据指针
 * @return 指向基本变量数据的指针
 */
Data* Equation::getPrim()
{
    return prim;
}

/**
 * @brief 获取守恒变量数据指针
 * @return 指向守恒变量数据的指针
 */
Data* Equation::getCons()
{
    return cons;
}

/**
 * @brief 获取右手端(RHS)数据指针
 * @return 指向右手端数据的指针
 */
Data* Equation::getRhs()
{
    return rhs;
}