#pragma once
#include "macro.hpp"

// 数据管理类，用于存储和操作计算数据
class Data {
public:
    // 构造函数
    Data() {};
    Data(int, int);
    Data(const Data&);

    // 初始化函数
    void solInit(int, int);
    void init(int, int);
    
    // 数据设置函数
    void setValue(std::vector<real>);
    void setvarName(std::vector<std::string>&&);
    
    // 数据访问操作符
    real& operator()(int, int);
    real& operator[](int);

    // 数据操作函数
    void operator+=(std::vector<real>);
    void setZeros();
    
    // 数据属性获取函数
    int size();
    int getNVar();
    int getN();

    // 数据提取函数
    std::vector<real> getIVar(int);
    
    // 迭代器函数
    std::vector<real>::iterator begin();
    std::vector<real>::iterator random(int i);
    std::vector<real>::iterator end();

    // 范数计算函数
    real maxElement(int);
    real getL2(int ivar);
    real getLinf(int ivar);
    
    // 变量名称
    std::vector<std::string> varName;

private:
    // 数据存储
    std::vector<real> data;  // 实际数据存储
    int n = 200;             // 数据大小
    int nVar = 1;            // 变量数量
};

// 数据读取器类，用于访问Data类中的特定数据
class DataReader {
public:
    // 构造函数
    DataReader() { }
    DataReader(int n_, int i0_, int offset_, int idim_, Data* data_)
        : n(n_)
        , offset(offset_)
        , data(data_)
        , idim(idim_)
    {
        i0 = i0_;
    }

    // 数据访问函数
    std::vector<real>::iterator operator()(int i)
    {
        return data->random(i0 + offset * i);
    }
    
    // 属性获取函数
    int getN() { return n; }
    int getNVar() { return data->getNVar(); }
    int getOffset() { return offset * data->getNVar(); }
    int idim;

protected:
    // 数据读取相关参数
    int n, i0, offset;
    Data* data;
};