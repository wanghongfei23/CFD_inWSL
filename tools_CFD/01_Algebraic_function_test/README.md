# CFD 代数函数测试工具

这是一个用于测试计算流体力学(CFD)中插值格式精度的工具。该项目通过使用特定的代数函数 `g_k(x) = x² * exp(0.75*(x-1))` 来评估插值算法的精度和收敛性。

## 项目概述

本项目主要用于验证高阶插值格式在CFD计算中的精度表现。通过在一个已知解析解的函数上测试插值格式，我们可以评估其逼近能力和收敛阶数。

### 测试函数

我们使用的测试函数为：
```
g_k(x) = x² * exp(0.75*(x-1))
```

该函数在CFD领域常被用来测试插值格式的精度，因为它具有适度的非线性特征，可以很好地模拟流体计算中遇到的实际问题。

## 项目结构

```
.
├── main.cpp                 # 主程序入口，执行插值精度测试
├── src/
│   ├── include/             # 头文件目录
│   │   ├── algebraic_function.hpp  # 代数函数定义及相关操作
│   │   ├── interscheme.hpp         # 插值格式实现
│   │   └── error_analysis.hpp      # 误差分析工具
│   └── src/                 # 源文件目录
│       ├── algebraic_function.cpp  # 代数函数实现
│       ├── error_analysis.cpp      # 误差分析实现
│       └── interscheme.cpp         # 插值格式实现
├── tests/
│   └── test_interpolation.cpp      # 插值测试程序
├── README.md                # 项目说明文档
└── Makefile/                # 构建系统配置目录
```

## 核心组件

### 1. 代数函数模块 (`algebraic_function.hpp/cpp`)

实现了测试函数 `g_k(x) = x² * exp(0.75*(x-1))` 及其一阶和二阶导数的计算。

### 2. 插值格式模块 (`interscheme.hpp/cpp`)

实现了三种不同的插值格式：
1. **whf_TCNS_AS_myF203_NoS** - 原始的高精度自适应插值格式
2. **weno5_JSchen** - 经典的五阶WENO格式
3. **Nicest5** - NICEST五阶格式

### 3. 误差分析模块 (`error_analysis.hpp/cpp`)

提供了多种误差度量方法：
- L1范数误差
- L2范数误差
- L∞范数误差（最大误差）
- 收敛阶计算

## 编译和运行

### 编译

```bash
# 编译主程序
g++ -std=c++11 main.cpp src/src/algebraic_function.cpp src/src/error_analysis.cpp src/src/interscheme.cpp -o cfd_test

# 编译测试程序
g++ -std=c++11 tests/test_interpolation.cpp src/src/algebraic_function.cpp src/src/error_analysis.cpp src/src/interscheme.cpp -o test_interp
```

### 运行

```bash
# 运行主程序（默认使用whf_TCNS_AS_myF203_NoS格式）
./cfd_test

# 运行主程序并指定插值格式
./cfd_test weno5_JSchen     # 使用WENO5_JSchen格式
./cfd_test whf_TCNS_AS_myF203_NoS  # 使用whf_TCNS_AS_myF203_NoS格式（默认）
./cfd_test Nicest5          # 使用Nicest5格式

# 或使用简写形式
./cfd_test weno5
./cfd_test whf
./cfd_test nicest5

# 运行测试程序（测试所有格式）
./test_interp
```

## 测试原理

程序在区间 [0.1, 1.0] 上使用不同分辨率的网格（10, 20, 40, 80, 160个网格点）进行测试：

1. 在每个网格点上计算测试函数的精确值
2. 使用五点模板插值格式计算数值解
3. 计算数值解与精确解之间的误差（L1、L2和L∞范数）
4. 根据不同网格分辨率下的误差计算收敛阶

## 添加新的插值格式

要添加新的插值格式，请按照以下步骤操作：

1. 在 [interscheme.hpp](file:///home/archwanghongfei/Documents/GitHub/CFD_codes_tools/tools_CFD/01_Algebraic_function_test/src/include/interscheme.hpp) 中声明新的插值函数
2. 在 [interscheme.cpp](file:///home/archwanghongfei/Documents/GitHub/CFD_codes_tools/tools_CFD/01_Algebraic_function_test/src/src/interscheme.cpp) 中实现具体的插值算法
3. 在 [main.cpp](file:///home/archwanghongfei/Documents/GitHub/CFD_codes_tools/tools_CFD/01_Algebraic_function_test/main.cpp) 中添加命令行参数识别逻辑
4. 在 [test_interpolation.cpp](file:///home/archwanghongfei/Documents/GitHub/CFD_codes_tools/tools_CFD/01_Algebraic_function_test/tests/test_interpolation.cpp) 中添加对该格式的测试

这种方法确保了良好的可扩展性，同时保持了代码的简洁性。

## 输出示例

程序会输出类似以下格式的结果：

```
CFD Algebraic Function Test
Testing interpolation scheme accuracy for g_k(x) = x^2 * exp(0.75*(x-1))
Using scheme: whf_TCNS_AS_myF203_NoS

 Grid Size  Grid Spacing      L1 Error      L2 Error     Linf Error
----------------------------------------------------------------------
        10        0.090000  1.0852e-03  1.3540e-03  3.5012e-03
        20        0.045000  6.8123e-05  8.5120e-05  2.2105e-04
        40        0.022500  4.2647e-06  5.3321e-06  1.3867e-05
        80        0.011250  2.6675e-07  3.3356e-07  8.6734e-07
       160        0.005625  1.6678e-08  2.0856e-08  5.4231e-08

Convergence Orders:
----------------------------------------
       Error Type               Order
             L1 Norm                4.00
             L2 Norm                4.00
            Linf Norm               4.00
```

## 应用场景

此工具可用于：
1. 验证新的插值格式的精度和收敛性
2. 比较不同插值格式的性能
3. CFD求解器中插值算法的验证和调试
4. 教学演示插值理论的实际应用

## 依赖项

- C++11或更高版本编译器
- 标准数学库（`<cmath>`）

## 注意事项

1. 当前实现中边界点采用了简单的直接赋值处理，在实际应用中可能需要更复杂的边界处理策略。
2. 插值格式目前固定使用五点模板，在近边界区域可能需要特殊处理。