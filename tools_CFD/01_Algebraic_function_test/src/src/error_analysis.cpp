#include "../include/error_analysis.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

/**
 * @brief 计算L1范数误差
 * @param exact 精确解
 * @param numerical 数值解
 * @return L1范数误差
 */
double ErrorAnalysis::l1_norm(const std::vector<double>& exact, 
                             const std::vector<double>& numerical) {
    if (exact.size() != numerical.size()) {
        throw std::invalid_argument("Exact and numerical solution vectors must have the same size");
    }

    double sum = 0.0;
    for (size_t i = 0; i < exact.size(); ++i) {
        sum += std::abs(exact[i] - numerical[i]);
    }
    return sum / exact.size();
}

/**
 * @brief 计算L2范数误差
 * @param exact 精确解
 * @param numerical 数值解
 * @return L2范数误差
 */
double ErrorAnalysis::l2_norm(const std::vector<double>& exact, 
                             const std::vector<double>& numerical) {
    if (exact.size() != numerical.size()) {
        throw std::invalid_argument("Exact and numerical solution vectors must have the same size");
    }

    double sum = 0.0;
    for (size_t i = 0; i < exact.size(); ++i) {
        double diff = exact[i] - numerical[i];
        sum += diff * diff;
    }
    return std::sqrt(sum / exact.size());
}

/**
 * @brief 计算L∞范数误差（最大误差）
 * @param exact 精确解
 * @param numerical 数值解
 * @return L∞范数误差
 */
double ErrorAnalysis::linf_norm(const std::vector<double>& exact, 
                               const std::vector<double>& numerical) {
    if (exact.size() != numerical.size()) {
        throw std::invalid_argument("Exact and numerical solution vectors must have the same size");
    }

    double max_error = 0.0;
    for (size_t i = 0; i < exact.size(); ++i) {
        double error = std::abs(exact[i] - numerical[i]);
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}

/**
 * @brief 计算收敛阶
 * @param errors 不同网格下的误差列表
 * @param grid_sizes 对应的网格大小
 * @return 收敛阶
 */
double ErrorAnalysis::convergence_order(const std::vector<double>& errors, 
                                       const std::vector<double>& grid_sizes) {
    if (errors.size() != grid_sizes.size() || errors.size() < 2) {
        throw std::invalid_argument("Errors and grid_sizes must have the same size and at least 2 elements");
    }

    // 使用最后两个网格点计算收敛阶
    size_t n = errors.size();
    double error_ratio = errors[n-2] / errors[n-1];
    double grid_ratio = grid_sizes[n-1] / grid_sizes[n-2];
    
    return std::log(error_ratio) / std::log(grid_ratio);
}