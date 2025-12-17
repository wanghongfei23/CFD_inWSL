#include "header.H"

/**
 * @brief 网格变量结构体
 * 
 * 该结构体封装了计算流体力学中常用的网格相关变量，
 * 包括网格坐标、解变量、通量和右端项等，用于含时间推进的计算
 */
struct GridVariables {
    GhostArray x;      // 网格坐标
    GhostArray u;      // 解变量
    GhostArray u_prev; // 前一时间步的解
    GhostArray flux;   // 通量存储在半节点上 (i+1/2位置)
    GhostArray rhs;    // 右端项

    real_t dx;         // 空间步长
    real_t time;       // 当前时间

    /**
     * @brief GridVariables构造函数（含时间推进版本）
     * @param N 网格点数
     * @param ghost 虚拟点数
     * @param length 计算域长度
     */
    GridVariables(int N, int ghost, real_t length)
        : x(N, 2 * ghost)
        , // N个网格点
        u(N, 2 * ghost)
        , u_prev(N, 2 * ghost)
        , flux(N, ghost)
        , // N个区间有N+1个界面
        rhs(N, 0)
        , // 右端项与解变量同尺寸
        dx(length / N)
        , time(0.0Q)
    {
        // 初始化网格坐标 (可选)
        for (int i : x.domain()) {
            x[i] = i * dx;
        }
        x.fill_periodic();
    }
};

/**
 * @brief 将计算结果写入Tecplot可视化文件
 * @param grid 网格变量结构体
 * @param current_time 当前时间
 * @param filename 输出文件名
 * 
 * 该函数将计算结果以Tecplot格式写入文件，便于后续可视化分析
 */
void write_to_tecplot_file(const GridVariables& grid, real_t current_time, const std::string& filename)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        throw std::runtime_error("Failed to open file for writing.");
    }

    outfile << "TITLE = \"Simulation Data\"\n";
    outfile << "VARIABLES = \"X\", \"U\"\n";
    outfile << "ZONE T=\"t=" << current_time << "\", I=" << grid.x.domain().length << "\n";

    for (int i : grid.x.domain()) {
        outfile << grid.x[i] << " " << grid.u[i] << "\n";
    }

    outfile.close();
}

/**
 * @brief 主函数（含时间推进）
 * 
 * 该程序用于计算线性对流方程的完整时间演化过程，
 * 使用三阶Runge-Kutta时间推进方法，并分析不同网格密度下的收敛性
 */
int main()
{
    // 设置OpenMP线程数
    omp_set_num_threads(10);
    // 网格密度序列（用于收敛性分析）
    std::vector<int> grid_densities = { 50, 100, 150, 200, 300, 400, 500 };
    // std::vector<int> grid_densities = { 50, 100};

    // 参数设置
    const int ghost = 5;     // 每侧虚拟点数
    const real_t L = 1.0Q;   // 计算域长度
    const real_t T = 1.0Q;   // 总模拟时间
    const real_t CFL = 0.01Q; // CFL数

    // 【王鸿飞】
    // 计算每个网格密度对应的网格步长
    std::vector<real_t> Lpergrid_densities;
    for (int density : grid_densities) 
    {
        Lpergrid_densities.push_back(L / density);
    }
    
    // 存储每种网格密度下的误差向量
    std::vector<real_t> L1_errors;
    std::vector<real_t> L2_errors;
    std::vector<real_t> Linf_errors;

    // 创建结果文件
    std::ofstream results_file("results.txt");
    if (!results_file.is_open()) {
        std::cerr << "Failed to open results.txt for writing." << std::endl;
        return 1;
    }

    // 写入文件头部
    results_file << "Variables = \"dx\", \"L1_Error\", \"L1_Order\", \"L2_Error\", \"L2_Order\", \"Linf_Error\", \"Linf_Order\"\n";
    results_file << "ZONE T=\"Convergence Analysis\"\n";

    // results_file << "N\tL1_Error\tL2_Error\tLinf_Error\tOrder_of_Accuracy\n";

    // 遍历所有网格密度进行计算
    for (int N : grid_densities) {
        // 创建网格变量集合
        GridVariables grid(N, ghost, L);

        // 初始化网格
// #pragma omp parallel for
        for (GhostArray::Range::Iterator i = grid.x.domain().begin(); i != grid.x.domain().end(); ++i) {
            grid.x[*i] = (*i + 0.5Q) * grid.dx; // 单元中心型网格坐标
            // 方波初始条件
            grid.u[*i] = exact_initial_solution(grid.x[*i], 0.0);
        }
        // 填充周期边界条件
        grid.x.fill_periodic();
        grid.u.fill_periodic();

        // 时间推进循环
        const real_t dt = CFL * grid.dx; // 基本时间步长
        while (grid.time < T) {
            // 计算自适应时间步长
            real_t current_dt = dt;
            if (grid.time + dt > T) {
                current_dt = T - grid.time; // 最后一步调整
            }

            // 保存前一时间步的解
// #pragma omp parallel for
            for (GhostArray::Range::Iterator i = grid.u.domain().begin(); i != grid.u.domain().end(); ++i) {
                grid.u_prev[*i] = grid.u[*i];
            }

            // 演化解到下一时间步（使用三阶Runge-Kutta方法）
            evolve(grid.u, grid.u_prev, grid.flux, grid.rhs, current_dt, grid.dx);

            // 更新当前时间
            grid.time += current_dt;
        }
        
        // 将结果写入Tecplot文件
        std::string filename = "output_N" + std::to_string(N) + ".dat";
        write_to_tecplot_file(grid, grid.time, filename);

        // 计算各种范数误差
        real_t L1_error = L1_norm_error(grid.u, grid.x, 0.0Q);
        real_t L2_error = L2_norm_error(grid.u, grid.x, 0.0Q);
        real_t Linf_error = Linf_norm_error(grid.u, grid.x, 0.0Q);

        // 存储误差
        L1_errors.push_back(L1_error);
        L2_errors.push_back(L2_error);
        Linf_errors.push_back(Linf_error);

        // 输出误差信息
        std::cout << "Grid density: " << N << "\n";
        std::cout << "L1 norm error: " << L1_error << "\n";
        std::cout << "L2 norm error: " << L2_error << "\n";
        std::cout << "L-infinity norm error: " << Linf_error << "\n\n";
    }

    // 计算收敛阶
    std::cout << "Order of accuracy:\n";
    // results_file << grid_densities[0] << "\t" << std::setprecision(5) << std::scientific 
    results_file << Lpergrid_densities[0] << "\t" << std::setprecision(5) << std::scientific 
                 << L1_errors[0] << "\t0\t" << L2_errors[0] << "\t0\t" << Linf_errors[0] << "\t0\n";
    for (size_t i = 1; i < grid_densities.size(); ++i) {
        real_t h1 = 1.0Q / grid_densities[i - 1];
        real_t h2 = 1.0Q / grid_densities[i];

        // 使用公式计算收敛阶：p = log(error1/error2) / log(h1/h2)
        real_t L1_order = log2(L1_errors[i - 1] / L1_errors[i]) / log2(h1 / h2);
        real_t L2_order = log2(L2_errors[i - 1] / L2_errors[i]) / log2(h1 / h2);
        real_t Linf_order = log2(Linf_errors[i - 1] / Linf_errors[i]) / log2(h1 / h2);

        std::cout << "From N = " << grid_densities[i - 1] << " to N = " << grid_densities[i] << ":\n";
        std::cout << "L1 norm order: " << L1_order << "\n";
        std::cout << "L2 norm order: " << L2_order << "\n";
        std::cout << "L-infinity norm order: " << Linf_order << "\n\n";

        // 更新收敛阶到结果文件
        // results_file << grid_densities[i] << "\t" << std::setprecision(5) << std::scientific
        results_file << Lpergrid_densities[i] << "\t" << std::setprecision(5) << std::scientific
                     << L1_errors[i] << "\t" << L1_order << "\t"
                     << L2_errors[i] << "\t" << L2_order << "\t"
                     << Linf_errors[i] << "\t" << Linf_order << "\n";
    }
    results_file.close();

    std::cout << "【with_time】: " << std::endl;
    std::cout << "使用的插值格式为: " << global_string << std::endl;
    
    return 0;
}