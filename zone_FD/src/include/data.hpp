#pragma once
#include "macro.hpp"

/**
 * @brief 数据管理类，用于存储和操作计算数据
 * 
 * Data类是一个通用的数据容器，用于存储CFD计算中的各种数据，
 * 如守恒变量、通量、源项等。它提供了一种类似矩阵的访问方式，
 * 支持按网格点和变量索引访问数据。
 */
class Data {
public:
    /**
     * @brief 默认构造函数
     */
    Data() {};
    
    /**
     * @brief 带参数的构造函数
     * @param n_ 网格点数量
     * @param nVar_ 变量数量
     */
    Data(int n_, int nVar_);
    
    /**
     * @brief 拷贝构造函数
     * @param origin_ 被拷贝的Data对象
     */
    Data(const Data& origin_);

    /**
     * @brief 解初始化函数，用于设置解的初始条件
     * @param n_ 网格点数量
     * @param nvar_ 变量数量
     */
    void solInit(int n_, int nvar_);
    
    /**
     * @brief 初始化函数
     * @param n_ 网格点数量
     * @param nvar_ 变量数量
     */
    void init(int n_, int nvar_);
    
    /**
     * @brief 设置数据值
     * @param value 包含所有数据的向量
     */
    void setValue(std::vector<real> value);
    
    /**
     * @brief 设置变量名称
     * @param varName_ 变量名称的向量（右值引用）
     */
    void setvarName(std::vector<std::string>&& varName_);
    
    /**
     * @brief 重载括号操作符，用于访问特定网格点和变量的数据
     * @param i 网格点索引
     * @param ivar 变量索引
     * @return 对应位置数据的引用
     */
    real& operator()(int i, int ivar);
    
    /**
     * @brief 重载方括号操作符，用于按线性索引访问数据
     * @param i 线性索引
     * @return 对应位置数据的引用
     */
    real& operator[](int i);

    /**
     * @brief 重载加等于操作符，将向量数据加到当前数据上
     * @param arr 要加上的数据向量
     */
    void operator+=(std::vector<real> arr);
    
    /**
     * @brief 将所有数据设置为零
     */
    void setZeros();
    
    /**
     * @brief 获取数据总大小
     * @return 数据总大小（网格点数×变量数）
     */
    int size();
    
    /**
     * @brief 获取变量数量
     * @return 变量数量
     */
    int getNVar();
    
    /**
     * @brief 获取网格点数量
     * @return 网格点数量
     */
    int getN();

    /**
     * @brief 提取特定变量的所有网格点数据
     * @param ivar 变量索引
     * @return 包含指定变量在所有网格点上数据的向量
     */
    std::vector<real> getIVar(int ivar);
    
    /**
     * @brief 获取数据起始迭代器
     * @return 指向第一个数据元素的迭代器
     */
    std::vector<real>::iterator begin();
    
    /**
     * @brief 获取指定网格点数据的迭代器
     * @param i 网格点索引
     * @return 指向第i个网格点数据的迭代器
     */
    std::vector<real>::iterator random(int i);
    
    /**
     * @brief 获取数据结束迭代器
     * @return 指向最后一个数据元素之后的迭代器
     */
    std::vector<real>::iterator end();

    /**
     * @brief 计算指定变量的最大值
     * @param ivar 变量索引
     * @return 指定变量的最大值
     */
    real maxElement(int ivar);
    
    /**
     * @brief 计算指定变量的L2范数
     * @param ivar 变量索引
     * @return 指定变量的L2范数
     */
    real getL2(int ivar);
    
    /**
     * @brief 计算指定变量的L∞范数
     * @param ivar 变量索引
     * @return 指定变量的L∞范数
     */
    real getLinf(int ivar);
    
    std::vector<std::string> varName; ///< 变量名称

private:
    std::vector<real> data;  ///< 实际数据存储
    int n = 200;             ///< 数据大小
    int nVar = 1;            ///< 变量数量
};

/**
 * @brief 数据读取器类，用于访问Data类中的特定数据
 * 
 * DataReader类提供了一个接口，用于从Data对象中读取特定部分的数据，
 * 主要用于处理多维数据或特定方向的数据切片。
 */
class DataReader {
public:
    /**
     * @brief 默认构造函数
     */
    DataReader() { }
    
    /**
     * @brief 带参数的构造函数
     * @param n_ 网格点数量
     * @param i0_ 起始索引
     * @param offset_ 偏移量
     * @param idim_ 维度标识
     * @param data_ 指向Data对象的指针
     */
    DataReader(int n_, int i0_, int offset_, int idim_, Data* data_)
        : n(n_)
        , offset(offset_)
        , data(data_)
        , idim(idim_)
    {
        i0 = i0_;
    }

    /**
     * @brief 重载函数调用操作符，用于获取指定网格点的迭代器
     * @param i 网格点索引
     * @return 指向第i个网格点数据的迭代器
     */
    std::vector<real>::iterator operator()(int i)
    {
        return data->random(i0 + offset * i);
    }
    
    /**
     * @brief 获取网格点数量
     * @return 网格点数量
     */
    int getN() { return n; }
    
    /**
     * @brief 获取变量数量
     * @return 变量数量
     */
    int getNVar() { return data->getNVar(); }
    
    /**
     * @brief 获取偏移量
     * @return 偏移量乘以变量数量
     */
    int getOffset() { return offset * data->getNVar(); }
    
    int idim; ///< 维度标识

protected:
    int n;      ///< 网格点数量
    int i0;     ///< 起始索引
    int offset; ///< 偏移量
    Data* data; ///< 指向Data对象的指针
};