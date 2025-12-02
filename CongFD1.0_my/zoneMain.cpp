/**
 * @file zoneMain.cpp
 * @brief 主程序入口文件，用于配置和运行CFD求解器
 */

#include "blockSolver.hpp"
#include "eigenSystem.hpp"
#include "000_globals.hpp"
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <limits>

#include <map>
#include <string>

/**
 * @brief 一维算例标识到字符串的映射表
 */
static std::map<int,std::string> exampleStr1D={
    {0,"Sod"},
    {1,"ShuOsher"},
    {2,"Lax"},
    {3,"sedov"},
    {4,"Woodward_Colella"},
    {5,"Double_sparse_wave"}
};

/**
 * @brief 二维算例标识到字符串的映射表
 */
static std::map<int,std::string> exampleStr2D={
    {0,"2D_Riemann_1"},
    {1,"2D_Riemann_2"},
    {2,"implosion"},
    {3,"RTI"},
    {4,"Double_Mach"},
    {5,"2D_Riemann_3"},
    {6,"KHI"}
};

/**
 * @brief 插值方法到字符串的映射表
 */
static std::map<InterMethod,std::string> disStr={
    {FIRSTORDER,"FIRSTORDER"},
    {MUSCL,"MUSCL"},
    {WCNS5,"WENO-JS"},
    {WCNSZ5,"WENO-Z"},
    {TCNS5,"TENO"},
    {WCNS5CONGZ,"TENO-S"},
    {WHFTCNSA,"TENO-myA"},
    {WHFTCNSASF002,"TENO-AS-myF002"},
    {WHFTCNSAH002,"TENO-AS-myH002"},
    {WHFTCNSASF102,"TENO-AS-myF102"},
    {WHFTCNSASF103,"TENO-AS-myF103"},
    {WHFTCNSASF102_reciprocal,"TENO-AS-myF102_reciprocal"},
    {WHFTCNSASF103_reciprocal,"TENO-AS-myF103_reciprocal"},
    {WHFTCNSAS_fx,"TENO-AS-fx"},
    {WHFTCNSAS_initial,"TENO-AS-initial"},
    {WHFTCNSAS_approx_1,"TENO-AS-approx_1"},
    {WHFTCNSAS_fx_real,"TENO-AS-fx-real"},
    {WHFTCNSAS_approx_2,"TENO-AS-approx_2"},
    {WHFTCNSASF202_2S,"TENO-AS-myF202_2S"},
    {WHFTCNSASF202_NoS,"TENO-AS-myF202_NoS"},
    {WHFTCNSASF203_NoS,"TENO-AS-myF203_NoS"},
    {temp011,"temp_name_011"},
    {temp012,"temp_name_012"},
    {temp013,"temp_name_013"},
    {temp014,"temp_name_014"},
    {temp015,"temp_name_015"},
    {temp016,"temp_name_016"},
    {temp017,"temp_name_017"},
    {temp018,"temp_name_018"},
    {temp019,"temp_name_019"}
};

/**
 * @brief 算例选择菜单
 */
void displayMenu() {
    std::cout << "\n=== CFD Solver Case Selection ===\n";
    std::cout << "1D Cases:\n";
    std::cout << "  0 - Sod\n";
    std::cout << "  1 - ShuOsher\n";
    std::cout << "  2 - Lax\n";
    std::cout << "  3 - Sedov\n";
    std::cout << "  4 - Woodward_Colella\n";
    std::cout << "  5 - Double_sparse_wave\n\n";
    
    std::cout << "2D Cases:\n";
    std::cout << "  10 - 2D_Riemann_1\n";
    std::cout << "  11 - 2D_Riemann_2\n";
    std::cout << "  12 - implosion\n";
    std::cout << "  13 - RTI\n";
    std::cout << "  14 - Double_Mach\n";
    std::cout << "  15 - 2D_Riemann_3\n";
    std::cout << "  16 - KHI\n\n";
    
    std::cout << "Enter your choice (0-5 for 1D, 10-16 for 2D): ";
}

/**
 * @brief 获取用户选择的算例
 * @return 用户选择的算例编号
 */
int getUserChoice() {
    int choice;
    while (!(std::cin >> choice) || 
           (choice < 0 || (choice > 5 && choice < 10) || choice > 16)) {
        std::cout << "Invalid input. Please enter a valid choice (0-5 for 1D, 10-16 for 2D): ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return choice;
}

/**
 * @brief 根据用户选择配置算例参数
 * @param[in,out] info Info对象指针
 * @param[in] choice 用户选择的算例编号
 */
void configureCase(Info* info, int choice) {
    switch(choice) {
        case 0: // Sod tube
            info->endStep = 20;
            info->outputDt = 0.01;
            info->CFL = 0.5;
            info->nCase = 0;
            info->calZone = { -0.5, 0.5, 0, 0, 0, 0 };
            info->iMax = { 201, 2, 2 };
            info->dim = 1;
            std::cout << "Configured: Sod tube problem\n";
            break;
            
        case 1: // Shu-Osher
            info->endStep = 1;
            info->CFL = 0.5;
            info->outputDt = 1.8;
            info->nCase = 1;
            info->calZone = {0, 10.0, 0, 0, 0, 0};
            info->iMax = {201, 2, 2};
            info->dim = 1;
            std::cout << "Configured: Shu-Osher problem\n";
            break;
            
        case 2: // Lax
            info->endStep = 14;
            info->outputDt = 0.01;
            info->CFL = 0.5;
            info->nCase = 2;
            info->calZone = { -0.5, 0.5, 0, 0, 0, 0 };
            info->iMax = { 201, 2, 2 };
            info->dim = 1;
            std::cout << "Configured: Lax problem\n";
            break;
            
        case 3: // Sedov
            info->endStep = 1;
            info->outputDt = 0.001;
            info->CFL = 0.5;
            info->nCase = 3;
            info->calZone = { -2, 2, 0, 0, 0, 0 };
            info->iMax = { 400, 2, 2 };
            info->dim = 1;
            std::cout << "Configured: Sedov problem\n";
            break;
            
        case 4: // Woodward-Colella
            info->endStep = 1;
            info->outputDt = 0.038;
            info->CFL = 0.1;
            info->nCase = 4;
            info->calZone = { 0, 1, 0, 0, 0, 0 };
            info->iMax = { 401, 2, 2 };
            info->dim = 1;
            std::cout << "Configured: Woodward-Colella problem\n";
            break;
            
        case 5: // Double sparse wave
            info->endStep = 100;
            info->outputDt = 0.01;
            info->CFL = 0.5;
            info->nCase = 5;
            info->calZone = { -5, 5, 0, 0, 0, 0 };
            info->iMax = { 401, 2, 2 };
            info->dim = 1;
            std::cout << "Configured: Double sparse wave problem\n";
            break;
            
        case 10: // 2D Riemann 1
            info->endStep = 1;
            info->outputDt = 0.8;
            info->CFL = 0.5;
            info->nCase = 0;
            info->calZone = { -0.5, 0.5, -0.5, 0.5, 0, 0 };
            info->iMax = { 401, 401, 2 };
            info->dim = 2;
            std::cout << "Configured: 2D Riemann 1 problem\n";
            break;
            
        case 11: // 2D Riemann 2
            info->endStep = 1;
            info->outputDt = 0.3;
            info->CFL = 0.5;
            info->nCase = 1;
            info->calZone = { -0.5, 0.5, -0.5, 0.5, 0, 0 };
            info->iMax = { 801, 801, 2 };
            info->dim = 2;
            std::cout << "Configured: 2D Riemann 2 problem\n";
            break;
            
        case 12: // Implosion
            info->endStep = 25;
            info->outputDt = 0.1;
            info->CFL = 0.5;
            info->nCase = 2;
            info->calZone = { -0.3, 0.3, -0.3, 0.3, 0, 0 };
            info->iMax = { 401, 401, 2 };
            info->dim = 2;
            std::cout << "Configured: Implosion problem\n";
            break;
            
        case 13: // RTI
            info->endStep = 1;
            info->outputDt = 1.95;
            info->CFL = 0.5;
            info->nCase = 3;
            info->calZone = { 0, 0.25, 0, 1, 0, 0 };
            info->iMax = { 101, 401, 2 };
            info->dim = 2;
            info->sourceType = GRAVITY;
            std::cout << "Configured: RTI problem\n";
            break;
            
        case 14: // Double Mach
            info->endStep = 20;
            info->outputDt = 0.01;
            info->CFL = 0.5;
            info->nCase = 4;
            info->calZone = { 0, 4, 0, 1, 0, 0 };
            info->iMax = { 801, 201, 2 };
            info->dim = 2;
            std::cout << "Configured: Double Mach problem\n";
            break;
            
        case 15: // 2D Riemann 3
            info->endStep = 1;
            info->outputDt = 0.25;
            info->CFL = 0.5;
            info->nCase = 5;
            info->calZone = { -0.5, 0.5, -0.5, 0.5, 0, 0 };
            info->iMax = { 801, 801, 2 };
            info->dim = 2;
            std::cout << "Configured: 2D Riemann 3 problem\n";
            break;
            
        case 16: // KHI
            info->endStep = 1;
            info->outputDt = 0.25; // Adjust as needed
            info->CFL = 0.5;
            info->nCase = 6;
            info->calZone = { -0.5, 0.5, -0.5, 0.5, 0, 0 }; // Adjust as needed
            info->iMax = { 401, 401, 2 }; // Adjust as needed
            info->dim = 2;
            std::cout << "Configured: KHI problem\n";
            break;
    }
}

/**
 * @brief 主函数，程序入口点
 * @return 程序退出状态码
 */
int main()
{
    omp_set_num_threads(10);
    // omp_set_num_threads(1);

    Info* info = new Info;

    info->eqType = EULER;
    info->spMethod = WCNS5;
    // 差分方法选项
        info->diffMethod = MND6;
        // info->diffMethod = TRAD6;

    // 插值方法选项
        // info->interMethod=LINEAR5;
        // info->interMethod=WCNSZ5Char;
        // info->BVD=true;
        // info->interMethod = NICEST5;
        // info->interMethod=WCNS5CONG;
        // info->sourceType=GRAVITY;

        // info->interMethod = WCNS5; //weno5_JSchen
        // info->interMethod = WCNSZ5; //weno5_Z
        info->interMethod = TCNS5; //Teno5_Z
        // info->interMethod = WCNS5CONGZ;//Teno5_CongZ
        // info->interMethod = WHFTCNSA;
        // info->interMethod = WHFTCNSASF002;
        // info->interMethod = WHFTCNSAH002;
        // info->interMethod = WHFTCNSASF102;
        // info->interMethod = WHFTCNSASF103;
        // info->interMethod = WHFTCNSASF102_reciprocal;
        // info->interMethod = WHFTCNSASF103_reciprocal;
        // info->interMethod = WHFTCNSAS_fx; // 不进行函数拟合，直接用原来近似的指数形式CT'
        // info->interMethod = WHFTCNSAS_initial; // 算CT'，原封不动的叠加A和S
        // info->interMethod = WHFTCNSAS_approx_1; // 算CT'，叠加A和S，分母近似掉-1
        // info->interMethod = WHFTCNSAS_fx_real; // 算
        // info->interMethod = WHFTCNSAS_approx_2; // 具体代入，叠加A和S，分母近似掉-1

        // info->interMethod = WHFTCNSASF202_2S;
        // info->interMethod = WHFTCNSASF202_NoS;
        // info->interMethod = WHFTCNSASF203_NoS;
        // info->interMethod = temp011;
        // info->interMethod = temp012;
        // info->interMethod = temp013;
        // info->interMethod = temp014;
        // info->interMethod = temp015;
        // info->interMethod = temp016;
        // info->interMethod = temp017;
        // info->interMethod = temp018;
        // info->interMethod = temp019;
        // info->interMethod = TCNSCongA;


    // 【预设算例编号】在此处设置需要运行的算例编号
    // 1D Cases: 
        // 0 - Sod, 
        // 1 - ShuOsher, 
        // 2 - Lax, 
        // 3 - Sedov, 
        // 4 - Woodward_Colella, 
        // 5 - Double_sparse_wave
    // 2D Cases: 
        // 10 - 2D_Riemann_1, 
        // 11 - 2D_Riemann_2, 
        // 12 - implosion, 
        // 13 - RTI, 
        // 14 - Double_Mach, 
        // 15 - 2D_Riemann_3, 
        // 16 - KHI

    const int presetCase = 2; // 算例选择
    
    // 移除交互式选择，直接使用预设的算例编号
    configureCase(info, presetCase);

    // 从info.txt文件中读取仿真参数配置
    std::ifstream file("info.txt");
    if (file.is_open()) {
        int n;
        real nf;
        file >> n;
        // 读取插值方法类型，并检查有效性
        if (n < temp019)
            info->interMethod = (InterMethod)n;

        file >> n;
        // 设置结束时间步
        info->endStep = n;

        file >> nf;
        // 设置输出时间间隔
        info->outputDt = nf;

        file >> nf;
        // 设置CFL数
        info->CFL = nf;

        file >> n;
        // 读取案例编号
        info->nCase = n;

        real nf1, nf2, nf3, nf4, nf5, nf6;
        file >> nf1;
        file >> nf2;
        file >> nf3;
        file >> nf4;
        file >> nf5;
        file >> nf6;
        // 设置计算区域范围
        info->calZone = { nf1, nf2, nf3, nf4, nf5, nf6 };

        int n1, n2, n3;
        file >> n1;
        file >> n2;
        file >> n3;
        // 设置网格点数
        info->iMax = { n1, n2, n3 };

        file >> n;
        // 设置问题维度
        info->dim = n;

        file >> n;
        // 设置OpenMP线程数
        omp_set_num_threads(n);

        std::cout << "file mode initialization finished\n";

    } else {
        // 文件打开失败提示
        std::cout << "file mode initialization failed\n";
    }

    InterMethod interscheme;

    // 创建BlockSolver对象
    BlockSolver bSolver(info);
    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();
    // 根据方程类型选择不同的时间步进循环方式
    if (info->eqType != EULER)
        bSolver.stepsLoop();
    else
        bSolver.stepsLoopCFL();
    // 结束计时并计算运行时间
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    // bSolver.stepsLoopDTS();
    // bSolver.solve();
    // bSolver.outputPrim();
    // bSolver.Test();

    // 命名系统
    // 查询命名映射表
    std::string caseName = "unknown";
    if (info->dim == 1) {
        auto it = exampleStr1D.find(info->nCase);
        if (it != exampleStr1D.end()) {
            caseName = it->second;
        }
    } else if (info->dim == 2) {
        auto it = exampleStr2D.find(info->nCase);
        if (it != exampleStr2D.end()) {
            caseName = it->second;
        }
    }
    // 添加网格信息到文件名
    std::string gridInfo = "";
    if (info->dim == 1) {
        gridInfo = std::to_string(info->iMax[0]);
    } else if (info->dim == 2) {
        gridInfo = std::to_string(info->iMax[0]) + "x" + std::to_string(info->iMax[1]);
    } else if (info->dim == 3) {
        gridInfo = std::to_string(info->iMax[0]) + "x" + std::to_string(info->iMax[1]) + "x" + std::to_string(info->iMax[2]);
    }
    
    std::string timeInfoFilename = caseName + " - " + disStr[info->interMethod] + " - " + gridInfo + ".txt";
    std::ofstream timeinfo(timeInfoFilename);
    
    // 时间信息打印与输出
    std::cout << "totaltime= " << duration << "   Finish\n";
    std::cout << "time= " << timepp / 1e6 << "   Finish\n";
    std::cout << "timesteps= " << bSolver.timesteps << "   Finish\n";
    std::cout << "solvertime= " << timesss << '\n';

    timeinfo << info->interMethod << std::endl;
    timeinfo << "totaltime= " << duration << "   Finish\n";
    timeinfo << "time= " << timepp / 1e6 << "   Finish\n";
    timeinfo << "timesteps= " << bSolver.timesteps << "   Finish\n";
    timeinfo << "solvertime= " << timesss << '\n';
    

    /// @brief 统计不同CT值的数量并输出到文件
    /// @details 当存在任何CT值计数时，将统计信息写入独立的文本文件
    if (global_counter_5 || global_counter_6 || global_counter_7 || global_counter_8 || global_counter_9 || global_counter_10) {
        std::string statFileName = caseName + " - " + disStr[info->interMethod] + " - " + gridInfo + "_CT-A_counter.txt";
        std::ofstream statFile(statFileName);
        statFile << "CT Value Statistics:\n";
        statFile << "CT=1e-05 count: " << global_counter_5 << "\n";
        statFile << "CT=1e-06 count: " << global_counter_6 << "\n";
        statFile << "CT=1e-07 count: " << global_counter_7 << "\n";
        statFile << "CT=1e-08 count: " << global_counter_8 << "\n";
        statFile << "CT=1e-09 count: " << global_counter_9 << "\n";
        statFile << "CT=1e-10 count: " << global_counter_10 << "\n";
        statFile.close();
    }
    /// @brief 检查是否存在非预期范围的CT值
    /// @details 如果pandaun_001为true，表示存在超出-5到-10范围的CT值
    if (pandaun_001)
    {
        std::cout << "有 非-5到-10的预期值" << "\n";
    }
    else
    {
        std::cout << "无 非-5到-10的预期值" << "\n";
    }
}