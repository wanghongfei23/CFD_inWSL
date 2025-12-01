#ifndef ERROR_ANALYSIS_HPP
#define ERROR_ANALYSIS_HPP

#include <vector>

/**
 * @brief 误差分析类，用于计算不同范数下的误差
 */
class ErrorAnalysis {
public:
    /**
     * @brief 计算L1范数误差
     * @param exact 精确解
     * @param numerical 数值解
     * @return L1范数误差
     */
    static double l1_norm(const std::vector<double>& exact, const std::vector<double>& numerical);

    /**
     * @brief 计算L2范数误差
     * @param exact 精确解
     * @param numerical 数值解
     * @return L2范数误差
     */
    static double l2_norm(const std::vector<double>& exact, const std::vector<double>& numerical);

    /**
     * @brief 计算L∞范数误差（最大误差）
     * @param exact 精确解
     * @param numerical 数值解
     * @return L∞范数误差
     */
    static double linf_norm(const std::vector<double>& exact, const std::vector<double>& numerical);

    /**
     * @brief 计算收敛阶
     * @param errors 不同网格下的误差列表
     * @param grid_sizes 对应的网格大小
     * @return 收敛阶
     */
    static double convergence_order(const std::vector<double>& errors, 
                                   const std::vector<double>& grid_sizes);
};

#endif // ERROR_ANALYSIS_HPP