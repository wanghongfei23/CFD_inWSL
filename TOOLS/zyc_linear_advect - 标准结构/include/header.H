// 包含所需的标准库和第三方库头文件
#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <quadmath.h>
#include <stdexcept>
#include <vector>

// 定义高精度浮点类型real_t为__float128
typedef __float128 real_t;

// 包含通量格式定义头文件
#include "flux_schemes.H"

/**
 * @brief 带虚拟点的数组类
 * 
 * 该类实现了一个带有虚拟点(Ghost Points)的数组结构，
 * 主要用于处理计算流体力学中的边界条件
 */
class GhostArray {
private:
    std::vector<real_t> data;  // 存储数据的向量
    int ghost_size;            // 虚拟点数量（每侧）
    int domain_size;           // 物理域内点数量

public:
    /**
     * @brief GhostArray构造函数
     * @param domain_points 物理域内的点数
     * @param ghost_points 每侧的虚拟点数
     */
    GhostArray(int domain_points, int ghost_points)
        : data(domain_points + 2 * ghost_points) // 必须先初始化data
        , ghost_size(ghost_points)               // 然后ghost_size
        , domain_size(domain_points)             // 最后domain_size
    {
    }

    /**
     * @brief GhostArray拷贝构造函数
     * @param other 被拷贝的对象
     */
    GhostArray(const GhostArray& other)
        : data(other.data)          // 保持与声明顺序一致
        , ghost_size(other.ghost_size)
        , domain_size(other.domain_size)
    {
    }

    /**
     * @brief 获取虚拟点数
     * @return 每侧的虚拟点数
     */
    int get_ghost_size() const { return ghost_size; }
    
    /**
     * @brief 获取物理域点数
     * @return 物理域内的点数
     */
    int get_domain_size() const { return domain_size; }

    /**
     * @brief 获取包括虚拟点在内的总大小
     * @return 总点数
     */
    int total_size() const { return data.size(); }

    /**
     * @brief 重载下标运算符（可修改版本）
     * @param i 索引（相对于物理域起始点）
     * @return 对应位置元素的引用
     */
    real_t& operator[](int i)
    {
        if (i < -ghost_size || i >= domain_size + ghost_size) {
            throw std::out_of_range("Index out of ghost array bounds");
        }
        return data[i + ghost_size];
    }

    /**
     * @brief 重载下标运算符（只读版本）
     * @param i 索引（相对于物理域起始点）
     * @return 对应位置元素的常量引用
     */
    const real_t& operator[](int i) const
    {
        if (i < -ghost_size || i >= domain_size + ghost_size) {
            throw std::out_of_range("Index out of ghost array bounds");
        }
        return data[i + ghost_size];
    }

    /**
     * @brief 获取原始数据指针（可修改版本）
     * @return 数据数组的指针
     */
    real_t* raw_data() { return data.data(); }
    
    /**
     * @brief 获取原始数据指针（只读版本）
     * @return 数据数组的常量指针
     */
    const real_t* raw_data() const { return data.data(); }

    /**
     * @brief 填充周期边界条件的虚拟点
     * 
     * 该函数将物理域边界的值复制到相对应的虚拟点中，
     * 实现周期性边界条件
     */
    void fill_periodic()
    {
        // 左侧虚拟点填充
        for (int i = -ghost_size; i < 0; ++i) {
            data[i + ghost_size] = data[domain_size + i + ghost_size];
        }
        // 右侧虚拟点填充
        for (int i = domain_size; i < domain_size + ghost_size; ++i) {
            data[i + ghost_size] = data[i - domain_size + ghost_size];
        }
    }

    /**
     * @brief 范围结构体
     * 
     * 用于表示GhostArray的不同访问范围，
     * 支持迭代器遍历
     */
    struct Range {
        /**
         * @brief 范围迭代器结构体
         */
        struct Iterator {
            int current;
            
            /**
             * @brief Iterator构造函数
             * @param pos 初始位置
             */
            Iterator(int pos)
                : current(pos)
            {
            }

            /**
             * @brief 解引用运算符
             * @return 当前位置
             */
            int operator*() const { return current; }
            
            /**
             * @brief 前置递增运算符
             * @return 自身引用
             */
            Iterator& operator++()
            {
                ++current;
                return *this;
            }
            
            /**
             * @brief 后置递增运算符
             * @return 递增前的副本
             */
            Iterator operator++(int)
            {
                Iterator temp = *this;
                ++(*this);
                return temp;
            }
            
            /**
             * @brief 前置递减运算符
             * @return 自身引用
             */
            Iterator& operator--()
            {
                --current;
                return *this;
            }
            
            /**
             * @brief 后置递减运算符
             * @return 递减前的副本
             */
            Iterator operator--(int)
            {
                Iterator temp = *this;
                --(*this);
                return temp;
            }
            
            /**
             * @brief 加法赋值运算符
             * @param n 增加的步数
             * @return 自身引用
             */
            Iterator& operator+=(int n)
            {
                current += n;
                return *this;
            }
            
            /**
             * @brief 加法运算符
             * @param n 增加的步数
             * @return 新的迭代器
             */
            Iterator operator+(int n) const
            {
                Iterator temp = *this;
                temp += n;
                return temp;
            }
            
            /**
             * @brief 减法赋值运算符
             * @param n 减少的步数
             * @return 自身引用
             */
            Iterator& operator-=(int n)
            {
                current -= n;
                return *this;
            }
            
            /**
             * @brief 减法运算符
             * @param n 减少的步数
             * @return 新的迭代器
             */
            Iterator operator-(int n) const
            {
                Iterator temp = *this;
                temp -= n;
                return temp;
            }
            
            /**
             * @brief 迭代器间距离运算符
             * @param other 另一个迭代器
             * @return 距离
             */
            int operator-(const Iterator& other) const
            {
                return current - other.current;
            }
            
            /**
             * @brief 等于比较运算符
             * @param other 另一个迭代器
             * @return 是否相等
             */
            bool operator==(const Iterator& other) const { return current == other.current; }
            
            /**
             * @brief 不等于比较运算符
             * @param other 另一个迭代器
             * @return 是否不等
             */
            bool operator!=(const Iterator& other) const { return !(*this == other); }
            
            /**
             * @brief 小于比较运算符
             * @param other 另一个迭代器
             * @return 是否小于
             */
            bool operator<(const Iterator& other) const { return current < other.current; }
            
            /**
             * @brief 大于比较运算符
             * @param other 另一个迭代器
             * @return 是否大于
             */
            bool operator>(const Iterator& other) const { return current > other.current; }
            
            /**
             * @brief 小于等于比较运算符
             * @param other 另一个迭代器
             * @return 是否小于等于
             */
            bool operator<=(const Iterator& other) const { return !(*this > other); }
            
            /**
             * @brief 大于等于比较运算符
             * @param other 另一个迭代器
             * @return 是否大于等于
             */
            bool operator>=(const Iterator& other) const { return !(*this < other); }
        };

        /**
         * @brief 获取起始迭代器
         * @return 起始迭代器
         */
        Iterator begin() const { return Iterator(start); }
        
        /**
         * @brief 获取结束迭代器
         * @return 结束迭代器
         */
        Iterator end() const { return Iterator(start + length); }

        int start;   // 起始索引
        int length;  // 范围长度
    };

    /**
     * @brief 获取物理域范围
     * @return 物理域范围对象
     */
    Range domain() const { return { 0, domain_size }; }
    
    /**
     * @brief 获取完整范围（包括虚拟点）
     * @return 完整范围对象
     */
    Range full_range() const { return { -ghost_size, domain_size + 2 * ghost_size }; }
};

/**
 * @brief 精确初始解函数
 * @param x 空间坐标
 * @param t 时间
 * @return 对应位置和时间的精确解
 * 
 * 该函数定义了线性对流方程的初始条件
 */

//  【王鸿飞】改函数
real_t exact_initial_solution(real_t x, real_t /*t*/)
{
    // Example: exact solution for a periodic sine function
    // u(x, t) = sin(2 * pi * (x - ct))
    // real_t c = 1.0Q; // wave speed
    // real_t x0 = x - c * t;
    // return sinq(2.0Q * M_PIq * x0) + 1.0Q; // Use M_PIq for quad-precision pi

    // Exact solution for a Gaussian function
    // u(x, t) = e^(-300 * (x - x_c)^2) where x_c = 0.5

    real_t x_c = 0.5Q;
    return expq(-300.0Q * (x - x_c) * (x - x_c));

    // return x * x * expq(0.75Q * (x - 1.0Q));

    // return x;
}

/**
 * @brief 精确导数函数
 * @param x 空间坐标
 * @param t 时间
 * @return 对应位置和时间的精确导数
 * 
 * 该函数定义了线性对流方程解的导数
 */
real_t exact_deriv(real_t x, real_t /*t*/)
{
    // real_t x_c = 0.5Q;
    // return 600.0Q * (x - x_c) * expq(-300.0Q * (x - x_c) * (x - x_c));
    // return -1.Q;
    return -expq(0.75Q * (x - 1.0Q)) * (2.0Q * x + 0.75Q * x * x);
}

/**
 * @brief 精确解函数
 * @param x 空间坐标
 * @param t 时间
 * @return 对应位置和时间的精确解
 * 
 * 该函数返回线性对流方程的精确解
 */
real_t exact_solution(real_t x, real_t /*t*/)
{
    // return exact_initial_solution(x, 0.0Q);
    return exact_deriv(x, 0.0Q);
}

/**
 * @brief 打印高精度实数
 * @param value 要打印的值
 * 
 * 该函数将__float128类型的数值转换为字符串并打印
 */
void print_real(real_t value)
{
    // Buffer for quadmath conversion
    char buf[128];

    // Convert quad precision to string with 10 decimal places
    quadmath_snprintf(buf, sizeof(buf), "%.10Qf", value);

    // Print to standard output
    std::cout << buf;
}

/**
 * @brief 重载输出流运算符
 * @param os 输出流
 * @param value 要输出的值
 * @return 输出流引用
 * 
 * 该函数允许直接将real_t类型的值输出到流中
 */
std::ostream& operator<<(std::ostream& os, const real_t& value)
{
    char buffer[128];
    quadmath_snprintf(buffer, sizeof(buffer), "%.30Qg", value);
    os << buffer;
    return os;
}

// 全局变量：是否使用周期边界条件
static bool is_periodic = true;

/**
 * @brief 计算右端项（WCNS5空间离散）
 * @param rhs 右端项数组
 * @param u 解变量数组
 * @param flux 通量数组
 * @param dx 空间步长
 * 
 * 该函数实现了WCNS5（Weighted Compact Nonlinear Scheme）格式的空间离散，
 * 计算线性对流方程的右端项
 */
void compute_rhs(GhostArray& rhs, GhostArray& u, GhostArray& flux, real_t dx)
{
    const real_t c = 1.0Q; // 波速
    if (is_periodic)
        u.fill_periodic();

    // 1. 一阶插值计算界面通量 (i+1/2位置)
// #pragma omp parallel for
    for (int i : flux.full_range()) {
        // 提取左值和右值
        // real_t u_left = u[i - 1];
        // real_t u_right = u[i];
        std::array<real_t, 5> q_left = { u[i - 3], u[i - 2], u[i - 1], u[i], u[i + 1] };
        std::array<real_t, 5> q_right = { u[i + 2], u[i + 1], u[i], u[i - 1], u[i - 2] };

        // 
        // 使用高阶插值得到左值和右值
        
        // 【王鸿飞】插值格式更换处
        // real_t u_left = whf_TcnsN_A(q_left);
        // real_t u_right = whf_TcnsN_A(q_right);
        // real_t u_left = whf_TcnsN_myASF203_NoS(q_left);
        // real_t u_right = whf_TcnsN_myASF203_NoS(q_right);
        real_t u_left = Linear_5th(q_left);
        real_t u_right = Linear_5th(q_right);

        // 一阶插值:取左右平均值
        // if (c > 0) {
        //     // c > 0, 使用前向差分
        //     flux[i] = c * u_left;
        // } else {
        //     // c <= 0, 使用后向差分
        //     flux[i] = c * u_right;
        // }

        real_t F_left = c * u_left;
        real_t F_right = c * u_right;
        flux[i] = (F_left + F_right) / 2.0Q - (c / 2.0Q) * (u_right - u_left);
    }
    if (is_periodic)
        flux.fill_periodic();
    std::array<real_t, 3> w = { 3.0Q / 2.0Q, -3.0Q / 10.0Q, 1.0Q / 30.0Q };
// std::array<real_t, 3> w = { 4.0Q / 3.0Q, -1.0Q / 6.0Q, 0.0Q };
// std::array<real_t, 3> w = { 1.0Q, 0.0Q, 0.0Q };
// 2. 二阶中心差分计算右端项
// #pragma omp parallel for
    for (int i : rhs.domain()) {
        // 二阶差分: (flux[i+1] - flux[i]) / dx
        // rhs[i] = -1.0Q / dx * (flux[i + 1] - flux[i]);
        rhs[i] = -w[0] / dx * (flux[i + 1] - flux[i])
            - w[1] / dx * (u[i + 1] - u[i - 1])
            - w[2] / dx * (flux[i + 2] - flux[i - 1]);
    }
}

/**
 * @brief 时间演化函数（三阶Runge-Kutta方法）
 * @param u 当前时间步的解
 * @param u_prev 前一时间步的解
 * @param flux 通量数组
 * @param rhs 右端项数组
 * @param dt 时间步长
 * @param dx 空间步长
 * 
 * 该函数使用三阶Runge-Kutta方法推进时间步，
 * 实现线性对流方程的时间积分
 */
void evolve(GhostArray& u, GhostArray& u_prev, GhostArray& flux, GhostArray& rhs,
    real_t dt, real_t dx) // 添加dx参数
{
    // 阶段1
    compute_rhs(rhs, u_prev, flux, dx); // 现在dx已定义
// #pragma omp parallel for
    for (int i : u.domain()) {
        u[i] = u_prev[i] + dt * rhs[i];
    }

    // 阶段2
    GhostArray temp(u);
    compute_rhs(rhs, u, flux, dx); // 传递dx
// #pragma omp parallel for
    for (int i : u.domain()) {
        temp[i] = 0.75Q * u_prev[i] + 0.25Q * (u[i] + dt * rhs[i]);
    }

    // 阶段3
    compute_rhs(rhs, temp, flux, dx); // 传递dx
// #pragma omp parallel for
    for (int i : u.domain()) {
        u[i] = (1.0Q / 3.0Q) * u_prev[i] + (2.0Q / 3.0Q) * (temp[i] + dt * rhs[i]);
    }
}

/**
 * @brief 计算L1范数误差
 * @param u 数值解
 * @param x 网格坐标
 * @param t 时间
 * @return L1范数误差
 * 
 * 该函数计算数值解与精确解之间的L1范数误差
 */
inline real_t L1_norm_error(const GhostArray& u, const GhostArray& x, real_t t)
{
    real_t error = 0.0Q;
    for (int i : u.domain()) {
        real_t xi = x[i];
        error += fabsq(u[i] - exact_solution(xi, t));
    }
    return error / u.get_domain_size();
}

/**
 * @brief 计算L2范数误差
 * @param u 数值解
 * @param x 网格坐标
 * @param t 时间
 * @return L2范数误差
 * 
 * 该函数计算数值解与精确解之间的L2范数误差
 */
inline real_t L2_norm_error(const GhostArray& u, const GhostArray& x, real_t t)
{
    real_t error = 0.0Q;
    for (int i : u.domain()) {
        real_t xi = x[i];
        error += powq(u[i] - exact_solution(xi, t), 2);
    }
    return sqrtq(error / u.get_domain_size());
}

/**
 * @brief 计算L∞范数误差
 * @param u 数值解
 * @param x 网格坐标
 * @param t 时间
 * @return L∞范数误差
 * 
 * 该函数计算数值解与精确解之间的L∞范数误差（最大误差）
 */
inline real_t Linf_norm_error(const GhostArray& u, const GhostArray& x, real_t t)
{
    real_t max_error = 0.0Q;
    for (int i : u.domain()) {
        real_t xi = x[i];
        real_t error = fabsq(u[i] - exact_solution(xi, t));
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}
