#pragma once
#include"data.hpp"
#include <cgnslib.h>

/**
 * @brief 网格块类，用于存储和操作计算网格信息
 * 
 * Block类负责管理计算域中的网格信息，包括网格顶点坐标、单元中心坐标以及相关的几何信息。
 * 它还提供了将网格数据输出为CGNS格式的功能。
 */
class Block
{
    public:
    /**
     * @brief 输出网格为CGNS格式文件
     * 
     * 将当前网格信息写入CGNS格式文件，便于后续可视化或与其他软件交互
     */
    void outputCgns();
    
    /**
     * @brief 获取单元中心坐标
     * 
     * 通过索引获取指定单元中心的坐标值
     * @param ic 单元索引
     * @param idim 坐标维度索引（0:x, 1:y, 2:z）
     * @return 指定单元在指定维度上的坐标值
     */
    real operator()(int ic,int idim);
    
    /**
     * @brief 获取内部单元最大索引数组
     * @return 包含各维度内部单元最大索引的数组
     */
    std::array<int,3> getICMax();
    
    /**
     * @brief 获取网格最大索引数组
     * @return 包含各维度网格最大索引的数组
     */
    std::array<int,3> getIMax();
    
    /**
     * @brief 获取指定维度的最小网格间距
     * @param i 维度索引
     * @return 在指定维度上的最小网格间距
     */
    real getMinDh(int i);
    
    /**
     * @brief 获取网格维度
     * @return 网格维度（1、2或3维）
     */
    int getDim();
    
    /**
     * @brief 获取单元中心坐标数据
     * @param idim 坐标维度索引（0:x, 1:y, 2:z）
     * @return 包含所有单元在指定维度上坐标值的向量
     */
    std::vector<real> getCellCoor(int idim);
    
    /**
     * @brief 获取顶点坐标数据
     * @param idim 坐标维度索引（0:x, 1:y, 2:z）
     * @return 包含所有顶点在指定维度上坐标值的向量
     */
    std::vector<real> getVertexCoor(int idim);
    
    /**
     * @brief 获取单元间距数据
     * @param idim 维度索引
     * @return 包含各单元在指定维度上网格间距的向量
     */
    std::vector<real> getCellInterval(int idim);

    protected:
    /// @cond
    // 友元类，允许Initializer访问私有成员
    friend class Initializer;
    /// @endcond
    
    // 网格数据成员
    Data coorVer;      ///< 网格顶点坐标数据
    Data coorCel;      ///< 网格单元中心坐标数据
    Data intervalCel;  ///< 网格单元间距数据
    
    // 网格基本信息
    int nVer,nCel;                    ///< 顶点数和单元数
    bool inited=false;                ///< 网格是否已初始化标志
    int dim;                          ///< 网格维度
    std::array<int,3> iMax,icMax;     ///< iMax通量点，icMax求解点

};