/**
 * @file zoneMain.cpp
 * @brief 主程序入口文件，用于配置和运行CFD求解器
 */

#include "blockSolver.hpp"
#include "eigenSystem.hpp"
#include "000_statistics_CTA.hpp"
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
    {WHFTCNSASF213_NoS,"TENO-AS-myF213_NoS"},
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
// =============================================================================
//                                omp线程数设置                                =
// =============================================================================
    omp_set_num_threads(10);
    // omp_set_num_threads(1);

// =============================================================================
//                                 Info 初始化                                 =
// =============================================================================
    Info* info = new Info;
    // --------------------------- 方程类型选项 --------------------------- 
    info->eqType = EULER;
    // --------------------------- 空间离散方法选项 --------------------------- 
    info->spMethod = WCNS5;
    // --------------------------- 差分方法选项 --------------------------- 
    info->diffMethod = MND6;
    // info->diffMethod = TRAD6;
    // --------------------------- 源项类型选项 --------------------------- 
    // info->sourceType=GRAVITY;
    // --------------------------- 插值方法选项 --------------------------- 
    // info->interMethod=LINEAR5;
    // info->interMethod=WCNSZ5Char;
    // info->interMethod = NICEST5;
    // info->interMethod=WCNS5CONG;

    // info->interMethod = WCNS5; //weno5_JSchen
    // info->interMethod = WCNSZ5; //weno5_Z
    // info->interMethod = TCNS5; //Teno5_Z
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
    info->interMethod = WHFTCNSASF213_NoS;
    // info->interMethod = temp012;
    // info->interMethod = temp013;
    // info->interMethod = temp014;
    // info->interMethod = temp015;
    // info->interMethod = temp016;
    // info->interMethod = temp017;
    // info->interMethod = temp018;
    // info->interMethod = temp019;
    // info->interMethod = TCNSCongA;

    // --------------------------- 算例选项 --------------------------- 
    const int presetCase = 
        // 0;  // Sod
        // 1;  // ShuOsher
        // 2;  // Lax
        // 3;  // Sedov
        // 4;  // Woodward_Colella
        // 5;  // Double_sparse_wave

        // 10; // 2D_Riemann_1
        // 11; // 2D_Riemann_2
        // 12; // implosion
        // 13; // RTI
        // 14; // Double_Mach
        15; // 2D_Riemann_3
        // 16; // KHI
    
    // 移除交互式选择，直接使用预设的算例编号
    configureCase(info, presetCase);

    // 从info.txt文件中读取仿真参数配置
    std::ifstream file("info.txt");
    if (file.is_open()) {
        int n;
        real nf;
        file >> n;                                    // 01. 读取插值方法类型，并检查有效性
        // if (n < temp019)
            info->interMethod = (InterMethod)n;
        file >> n;                                    // 02. 设置结束时间步
            info->endStep = n;
        file >> nf;                                   // 03. 设置输出时间间隔
            info->outputDt = nf;
        file >> nf;                                   // 04. 设置CFL数
            info->CFL = nf;
        file >> n;                                    // 05. 读取案例编号
            info->nCase = n;
        real nf1, nf2, nf3, nf4, nf5, nf6;
        file >> nf1;                                  // 06. 设置计算区域范围
        file >> nf2;                                  // 07. 设置计算区域范围
        file >> nf3;                                  // 08. 设置计算区域范围
        file >> nf4;                                  // 09. 设置计算区域范围
        file >> nf5;                                  // 10. 设置计算区域范围
        file >> nf6;                                  // 11. 设置计算区域范围
            info->calZone = { nf1, nf2, nf3, nf4, nf5, nf6 };
        int n1, n2, n3;
        file >> n1;                                   // 12. 设置网格点数
        file >> n2;                                   // 13. 设置网格点数
        file >> n3;                                   // 14. 设置网格点数
            info->iMax = { n1, n2, n3 };
        file >> n;                                    // 15. 设置问题维度
            info->dim = n;
        file >> n;                                    // 16. 设置OpenMP线程数
        omp_set_num_threads(n);

        std::cout << "file mode initialization finished（文件模式初始化完成）\n";

    } else {
        // 文件打开失败提示
        std::cout << "file mode initialization failed（文件模式初始化失败）\n";
    }

    // 创建BlockSolver对象，名字为bSolver，参数为info指针
    BlockSolver bSolver(info);

    // 开始计时（开始时间变量名为start）
    auto start = std::chrono::high_resolution_clock::now();

    // 根据方程类型选择不同的时间步进循环方式
    if (info->eqType != EULER)
        bSolver.stepsLoop();
    else
        bSolver.stepsLoopCFL();

    // 结束计时并计算运行时间（结束时间变量名为stop，持续时间变量名为duration）
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

    // output_CTA();

}
