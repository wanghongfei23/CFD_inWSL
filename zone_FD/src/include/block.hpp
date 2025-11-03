#pragma once
#include"data.hpp"
#include <cgnslib.h>

// 网格块类，用于存储和操作计算网格信息
class Block
{
    public:
    // 网格输出相关函数
    void outputCgns();
    
    // 坐标获取函数
    real operator()(int,int);
    
    // 网格索引相关函数
    std::array<int,3> getICMax();
    std::array<int,3> getIMax();
    
    // 网格属性获取函数
    real getMinDh(int i);
    int getDim();
    
    // 坐标数据获取函数
    std::vector<real> getCellCoor(int idim);
    std::vector<real> getVertexCoor(int idim);
    std::vector<real> getCellInterval(int);

    protected:
    // 友元类，允许Initializer访问私有成员
    friend class Initializer;
    
    // 网格数据成员
    Data coorVer;      // 网格顶点坐标数据
    Data coorCel;      // 网格单元中心坐标数据
    Data intervalCel;  // 网格单元间距数据
    
    // 网格基本信息
    int nVer,nCel;           // 顶点数和单元数
    bool inited=false;       // 网格是否已初始化标志
    int dim;                 // 网格维度
    std::array<int,3> iMax,icMax;  // iMax通量点，icMax求解点

};
