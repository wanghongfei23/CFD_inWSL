// 预处理指令
#pragma once            // 只包含一次该头文件，防止重复包含
// 宏定义
#define real double     // 定义 real 为 double 类型，便于后续统一修改精度
#define ind int         // 定义 ind 为 int 类型，便于统一索引类型
// 引入头文件
#include <array>        // 提供 std::array 容器
#include <cstring>      // 提供 C 风格字符串操作函数
#include <iostream>     // 提供标准输入输出流
#include <memory>       // 提供智能指针相关功能
#include <string>       // 提供 std::string 字符串类
#include <vector>       // 提供 std::vector 动态数组容器

#include <algorithm>    // 提供常用算法，如排序、查找等
#include <cassert>      // 提供断言功能
#include <cmath>        // 提供数学函数，如三角、指数等
#include <format>       // 提供格式化输出（C++20）
#include <map>          // 提供 std::map 映射容器

#include <chrono>       // 提供时间和计时相关功能

// 全局变量
inline long timepp = 0;   // 记录数值计算核心部分的累计执行时间
inline long timesss = 0;  // 记录求解器的累计执行时间

// 【RTI改gamma】
//其他算例的gamma
#define GAMMA 1.4
//RT算例的gamma
// #define GAMMA 5.0/3.0

// 边界条件类型
enum BndType {TYPENULL,PERIODIC1D,SYMMETRY1D,DIRICLET,DIRICLET_SODL,DIRICLET_SODR,FLUXGHOST,SUPERSONICOUTLET,DoubleMachUp,ARDBND,
    SYMMETRYX, SYMMETRYY}; //2D
// 【王鸿飞】begin插值格式enum
// 插值格式
enum InterMethod {FIRSTORDER,
    MUSCL,
    WCNS5,
    WCNSZ5,
    WCNS5Char,
    WCNSZ5Char,
    WCNS5CONG,
    TCNSCongA,
    WCNS5CONGZ,
    WCNS5CONGZCT4,
    WCNS5CONGZCT7,
    TCNS5,
    TCNS5CT4,
    TCNS5CT7,
    LINEAR5,
    MUCSLIN5,
    INTERMAX,
    
    whfTCNSN,
    whfTCNSNA,
    whfTCNSNAS,
    whfTCNSNS,
    whfTCNSNLAD,
    congTCNS5CT5,
    congTCNS5CT10,
    whfAITCNSNS,
    whfAITCNSNA,
    whfAITCNSNAS_1,
    whfAITCNSNAS_2,
    whfAITCNSNLADS,
    whfAITCNSNAZS,
    whfAITCNSNmyASF002_1,
    whfAITCNSNmyASF002_2,
    whfAITCNSNmyASF002_ai1,
    whfzycTCNSNmyASF002_1,
    whfCOMPARE}; 
// 【王鸿飞】end

// 插分格式
enum DiffMethod { HDS6,TRAD6,TRAD2,MND6 };
// 方程类型
enum EquationType { LINEARCONV1D,BURGERS1D,EULER,ACCURACYTEST };
// 通量构造方法
enum FluxMethod { HLLC,ROE };
// 时间离散格式
enum TimeMethod { RK3SSP,EulerFront };
// 源项类型
enum SourceType { SOURCENULL,GRAVITY };

// inline 内联函数，将三维数组的下标映射到一维数组的线性索引中
inline int index(int i, int j, int k, std::array<int, 3> iMax)
{
    return i + j * iMax[0] + k * iMax[0] * iMax[1];
}

/*
用于边界数据访问，可以实现不同的边界格式，比如周期性或者对称

计算三维数据在指定方向上的索引偏移量
输入：方向idim (1-y,2-z,3-x); 坐标i,j; 各维数组大小iMax
输出：返回一个包含起始索引和偏移量的数组
*/
inline std::array<int, 2> calOffset(int idim, int i, int j, std::array<int, 3> iMax)
{
    std::array<int, 3> offsets { 1, iMax[0], iMax[0] * iMax[1] };
    std::array<int, 2> res;
    if (idim == 1) {
        res[0] = i * offsets[1] + j * offsets[2]; // i0
        res[1] = offsets[0]; // offset
    } else if (idim == 2) {
        res[0] = i * offsets[0] + j * offsets[2]; // i0
        res[1] = offsets[1]; // offset
    } else if (idim == 3) {
        res[0] = i * offsets[0] + j * offsets[1]; // i0
        res[1] = offsets[2]; // offset
    }
    return res;
}
// 计算三维数据在指定方向上的索引偏移量（反向）
inline std::array<int, 2> calOffsetInverse(int idim, int i, int j,std::array<int, 3> iMax)
{
    std::array<int, 3> offsets { 1, iMax[0], iMax[0] * iMax[1] };
    std::array<int, 2> res;
    if (idim == 1) {
        res[0] = i * offsets[1] + j * offsets[2] + (iMax[0] - 1) * offsets[0]; // i0
        res[1] = -offsets[0]; // offset
    } else if (idim == 2) {
        res[0] = i * offsets[0] + j * offsets[2] + (iMax[1] - 1) * offsets[1]; // i0
        res[1] = -offsets[1]; // offset
    } else if (idim == 3) {
        res[0] = i * offsets[0] + j * offsets[1] + (iMax[2] - 1) * offsets[2]; // i0
        res[1] = -offsets[2]; // offset
    }
    return res;
}