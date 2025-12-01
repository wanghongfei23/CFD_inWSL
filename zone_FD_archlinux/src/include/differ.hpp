#pragma once

#include "data.hpp"

#include "info.hpp"
#include <functional>

// 差分基类，用于计算空间导数
class Differ {
public:
    // 构造函数
    Differ() {};
    Differ(std::shared_ptr<Data> fluxesHalf_, std::shared_ptr<Data> fluxesNode_,
        DataReader rhsR_)
        : fluxesHalf(fluxesHalf_)
        , fluxesNode(fluxesNode_)
        , rhsR(rhsR_) {};
    
    // 初始化函数
    void init(std::shared_ptr<Data> fluxesHalf_,
        std::shared_ptr<Data> fluxesNode_, DataReader rhsR_)
    {
        fluxesHalf = fluxesHalf_;
        fluxesNode = fluxesNode_;
        rhsR = rhsR_;
    }
    
    // 设置常数函数
    void setConstantH(real constantH_) { constantH = constantH_; }
    
    // 检查函数
    void check() { std::cout << "initialized successfully Differ\n"; }
    
    // 求解接口函数
    virtual void solve();

protected:
    // 数据成员
    std::shared_ptr<Data> fluxesHalf;  // 半节点通量数据
    std::shared_ptr<Data> fluxesNode;  // 节点通量数据
    DataReader rhsR;                   // 右手端数据读取器
    real constantH;                    // 常数
    int idim;                          // 维度索引
};

// 显式六阶差分类
class ExplicitDif6 : public Differ {
public:
    // 求解函数
    void solve() final;
};

// 中点和节点六阶差分类
class MidNodeAndNodeDif6 : public Differ {
public:
    // 求解函数
    void solve() final;
};