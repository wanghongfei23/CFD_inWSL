/**
 * @file test_interpolation.cpp
 * @brief 测试插值格式对代数函数的逼近效果
 */

#include "../src/include/algebraic_function.hpp"
#include "../src/include/interscheme.hpp"
#include "../src/include/error_analysis.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>

// 定义函数指针类型
typedef double (*InterpFunc)(std::array<double, 5>);

void test_interpolation_accuracy(InterpFunc interp_func, const std::string& scheme_name) {
    std::cout << "=== 插值格式精度测试 (" << scheme_name << ") ===\n";
    std::cout << "测试函数: g_k(x) = x^2 * exp(0.75*(x-1))\n\n";

    // 定义计算区间
    const double x_start = 0.1;
    const double x_end = 1.0;
    
    // 定义不同的网格分辨率进行收敛性测试
    std::vector<int> grid_sizes = {10, 20, 40, 80, 160};
    
    // 存储各网格分辨率下的误差
    std::vector<double> l1_errors, l2_errors, linf_errors;
    std::vector<double> grid_spacings;
    
    std::cout << std::setw(10) << "Grid Size" 
              << std::setw(15) << "Grid Spacing" 
              << std::setw(15) << "L1 Error" 
              << std::setw(15) << "L2 Error" 
              << std::setw(15) << "Linf Error" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    // 对每种网格分辨率进行测试
    for (int grid_size : grid_sizes) {
        const double dx = (x_end - x_start) / grid_size;
        grid_spacings.push_back(dx);
        
        // 创建网格点
        std::vector<double> x_points(grid_size + 1);
        std::vector<double> exact_values(grid_size + 1);
        
        for (int i = 0; i <= grid_size; ++i) {
            x_points[i] = x_start + i * dx;
            exact_values[i] = AlgebraicFunction::g_k(x_points[i]);
        }
        
        // 使用插值格式计算数值解
        std::vector<double> numerical_values(grid_size + 1);
        
        // 内部点使用5点模板插值
        for (int i = 2; i < grid_size - 2; ++i) {
            std::array<double, 5> stencil = {
                exact_values[i-2], exact_values[i-1], exact_values[i], exact_values[i+1], exact_values[i+2]
            };
            
            // 应用插值格式计算该点的值
            numerical_values[i] = interp_func(stencil);
        }
        
        // 边界点采用简单处理（实际应用中可能需要更复杂的边界处理）
        numerical_values[0] = exact_values[0];
        numerical_values[1] = exact_values[1];
        numerical_values[grid_size-1] = exact_values[grid_size-1];
        numerical_values[grid_size] = exact_values[grid_size];
        
        // 计算误差
        double l1_error = ErrorAnalysis::l1_norm(exact_values, numerical_values);
        double l2_error = ErrorAnalysis::l2_norm(exact_values, numerical_values);
        double linf_error = ErrorAnalysis::linf_norm(exact_values, numerical_values);
        
        l1_errors.push_back(l1_error);
        l2_errors.push_back(l2_error);
        linf_errors.push_back(linf_error);
        
        std::cout << std::setw(10) << grid_size 
                  << std::setw(15) << std::fixed << std::setprecision(6) << dx
                  << std::setw(15) << std::scientific << std::setprecision(6) << l1_error
                  << std::setw(15) << l2_error
                  << std::setw(15) << linf_error << std::endl;
    }
    
    // 计算收敛阶
    std::cout << "\n收敛阶:\n";
    std::cout << std::string(40, '-') << std::endl;
    std::cout << std::setw(20) << "误差类型" << std::setw(20) << "收敛阶" << std::endl;
    
    double l1_order = ErrorAnalysis::convergence_order(l1_errors, grid_spacings);
    double l2_order = ErrorAnalysis::convergence_order(l2_errors, grid_spacings);
    double linf_order = ErrorAnalysis::convergence_order(linf_errors, grid_spacings);
    
    std::cout << std::setw(20) << "L1 范数" << std::setw(20) << std::fixed << std::setprecision(2) << l1_order << std::endl;
    std::cout << std::setw(20) << "L2 范数" << std::setw(20) << l2_order << std::endl;
    std::cout << std::setw(20) << "Linf 范数" << std::setw(20) << linf_order << std::endl;
    
    std::cout << std::endl;
}

int main() {
    // 测试所有插值格式
    test_interpolation_accuracy(weno5_JSchen, "weno5_JSchen");
    test_interpolation_accuracy(whf_TCNS_AS_myF203_NoS, "whf_TCNS_AS_myF203_NoS");
    test_interpolation_accuracy(Nicest5, "Nicest5");
    
    return 0;
}