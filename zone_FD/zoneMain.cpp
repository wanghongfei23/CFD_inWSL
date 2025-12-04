
#include "blockSolver.hpp"    // 自定义求解器核心
#include "eigenSystem.hpp"    // 特征系统计算
#include <fstream>            // 文件输入输出流

int main()
{
    // 特征系统变换测试代码
    // 用于验证eigensystemEuler2D类的特征值变换功能正确性
    // 测试流程：创建测试用的原始变量 -> 执行特征变换 -> 执行逆变换 -> 验证结果准确性检测
    
    // std::array<real,4> prim0={1.0,0.75,-0.5,1.0};
    // std::array<real,4> prim={2.0,-0.75,0.5,1.0};
    // eigensystemEuler2D eig=eigensystemEuler2D(prim0,{1,0,0});
    // auto eigValues=eig.primToChar(prim);
    // auto prim2=eig.charToPrim(eigValues);
    // std::cout<<"finish\n";

    // ===================================求解器配置===================================

    // omp_set_num_threads(18);  // 设置OpenMP线程数
    omp_set_num_threads(10);  // 设置OpenMP线程数
    Info* info = new Info;  // 动态分配配置对象

    info->eqType = EULER;
    info->spMethod = WCNS5;
    info->fluxMethod = HLLC;
    info->diffMethod = MND6;
    // info->diffMethod = TRAD6;
    // info->BVD=true;
    // info->sourceType=GRAVITY;   

    // ------------------------interMethod--------------------
    // info->interMethod=WCNS5CONGZCT7;
    // info->interMethod = TCNS5;
    // info->interMethod=WCNSZ5;
    // info->interMethod=WCNS5Char;ROE
    // info->interMethod=WCNS5CONG;
    // info->interMethod=TCNSCongA;
    // info->interMethod=WCNS5CONGZ;
    // info->interMethod = WCNS5;    // weno5_JSchen
    
    // 【王鸿飞】begin插值格式开发（1.0）手搓
    // info->interMethod = whfTCNSN; //【通过】【待定】与师兄原有的 congTCNS5CT5 效果相同，验证一致
    // info->interMethod = whfTCNSNS; //作废
    // info->interMethod = whfTCNSNAS;//作废
    // info->interMethod = whfTCNSNLAD;//作废
    // 【王鸿飞】end插值格式开发（1.0）手搓

    // 【王鸿飞】begin插值格式开发（2.0）AI-原格式
    // info->interMethod= congTCNS5CT5 ; //【通过】师兄自带的，暂定为正确
    // info->interMethod= congTCNS5CT10 ; //【通过】师兄自带的，暂定为正确
    // info->interMethod= whfAITCNSNA ;
    // info->interMethod= whfAITCNSNS ;
    // info->interMethod= whfAITCNSNAS_1 ;
    // info->interMethod= whfAITCNSNAS_2 ;
    // info->interMethod= whfAITCNSNLADS ;
    // info->interMethod= whfAITCNSNAZS ;
    // 【王鸿飞】end插值格式开发（2.0）AI-原格式

    // 【王鸿飞】begin插值格式开发（3.0）AI-新格式
    // info->interMethod= whfAITCNSNmyASF002_1 ;
    // info->interMethod= whfAITCNSNmyASF002_2 ;
    // info->interMethod= whfAITCNSNmyASF002_ai1 ;
    // 【王鸿飞】end插值格式开发（3.0）AI-新格式

    // 【王鸿飞】begin插值格式（4.0）手搓及优化
    // info->interMethod= whfzycTCNSNmyASF002_1 ; //
    // info->interMethod = whfTCNSNA; //
    // info->interMethod = whfCOMPARE; //找不同
    // 【王鸿飞】end插值格式（4.0）手搓及优化

    // 【王鸿飞】begin插值格式（F203）)手搓及优化
    // info->interMethod= new_TCNS5 ;
    // info->interMethod= new_WHFTCNSA ;
    info->interMethod= new_WCNS5CONGZ ;
    // info->interMethod= new_WHFTCNSASF203_NoS ;
    // info->interMethod= new_WHFTCNSASF202_NoS ;

    // 【王鸿飞】end插值格式（F203）)手搓及优化

    // ===================================算例设置===================================

    // Shu-Osher
    // info->endStep=1;
    // info->CFL=0.5;
    // info->outputDt=1.8;
    // info->nCase=1;
    // info->calZone={0,10.0,0,0,0,0};
    // info->iMax={201,2,2};
    // info->dim=1;

    // sod tube
    // info->CFL = 0.5;
    // info->endStep = 2;
    // info->outputDt = 0.1;
    // info->nCase = 0;
    // info->calZone = { -0.5, 0.5, 0, 0, 0, 0 };
    // info->iMax = { 101, 2, 2 };
    // // info->iMax = { 1001, 2, 2 };
    // info->dim = 1;

    // lax sod tube
    // info->endStep=14;
    // info->outputDt=0.01;
    // info->CFL=0.5;
    // info->nCase=2;
    // info->calZone={-0.5,0.5,0,0,0,0};
    // info->iMax={201,2,2};
    // // info->iMax={2001,2,2};
    // info->dim=1;

    // lax sod tube speed test
    // info->endStep=14;
    // info->outputDt=0.01;
    // info->CFL=0.1;
    // info->nCase=2;
    // info->calZone={-0.5,0.5,0,0,0,0};
    // info->iMax={2001,2,2};
    // info->dim=1;

    // sedov
    // info->endStep=1;
    // info->outputDt=0.001;
    // info->CFL=0.5;
    // info->nCase=3;
    // info->calZone={-2,2,0,0,0,0};
    // info->iMax={400,2,2};
    // info->dim=1;

    // Woodward-Colella
    // info->endStep=1;
    // info->outputDt=0.038;
    // info->CFL=0.4;
    // // info->CFL=0.1;
    // info->nCase=4;
    // info->calZone={0,1,0,0,0,0};
    // // info->iMax={41,2,2};
    // info->iMax={401,2,2};
    // // info->iMax={4001,2,2};
    // info->dim=1;

    // 双稀疏波  Double sparse wave

    // info->endStep=100;
    // info->outputDt=0.01;

    // info->endStep=1;
    // info->outputDt=1.0;
    // info->CFL=0.5;
    // info->nCase=5;
    // info->calZone={-5,5,0,0,0,0};
    // info->iMax={401,2,2};
    // // info->iMax={11,2,2};
    // info->dim=1;

    // implosion
    //chenyuqing: 算例三：终极算例

    // info->endStep=1;
    // info->outputDt=0.1;
    // info->CFL=0.5;
    // info->nCase=2;
    // info->calZone={-0.3,0.3,-0.3,0.3,0,0};
    // info->iMax={401,401,2};
    // info->dim=2;

    // Riemann 1
    // /chenyuqing: 算例一

    // info->endStep = 1; //输出多少步
    // info->outputDt = 0.8; //步与步之间的间隔
    // info->CFL = 0.5;//CFL数
    // info->nCase = 0;
    // info->calZone = { -0.5, 0.5, -0.5, 0.5, 0, 0 };//计算域
    // info->iMax = { 401, 401, 2 };//网格数
    // info->dim = 2;


    // Riemann 2 vortex

    // info->endStep=1;
    // info->outputDt=0.3;
    // info->CFL=0.5;
    // info->nCase=1;
    // info->calZone={-0.5,0.5,-0.5,0.5,0,0};
    // info->iMax={801,801,2};
    // info->dim=2;

    // Riemann 3 
    //chenyuqing: 算例二
    
    // info->endStep = 1;
    // info->outputDt = 0.4;
    // info->CFL = 0.5;
    // info->nCase = 5;
    // info->calZone = { -0.5, 0.5, -0.5, 0.5, 0, 0 };
    // info->iMax = { 801, 801, 2 };
    // info->dim = 2;

    // RT instability
    // 记得改GAMMA (macro.hpp中)
    //chenyuqing: 算例四

    //  info->endStep=1;
    //  info->outputDt=1.95;
    //  info->CFL=0.5;
    //  info->nCase=3;
    //  info->calZone={0,0.25,0,1,0,0};
    // //  info->iMax={301,1201,2};
    // //  info->iMax={601,2401,2};
    //  info->iMax={65,257,2};
    //  info->dim=2;
    //  info->sourceType=GRAVITY;



    // info->diffMethod=HDS6;

    // Double Mach

    info->endStep=20;
    info->outputDt=0.01;
    
    // info->endStep=2;
    // info->outputDt=0.1;
    
    info->CFL=0.5;
    info->nCase=4;
    info->calZone={0,4,0,1,0,0};
    info->iMax={801,201,2};
    info->dim=2;


    // 【王鸿飞】begin新算例
    // KHI

    // info->endStep=10;
    // info->outputDt=0.1;
    // info->CFL=0.5;
    // // info->CFL=0.1;
    // info->nCase=6;
    // info->calZone={-0.5,0.5,-0.5,0.5,0,0};
    // info->iMax={257,257,2};
    // // info->iMax={513,513,2};
    // // info->iMax={1025,1025,2};
    // info->dim=2;

    // 【王鸿飞】end

    // 从info.txt配置文件中读取求解器运行参数
    // 配置参数优先级: 配置文件 > 硬编码默认值
    std::ifstream file("info.txt");
    if (file.is_open()) {
        int n;
        real nf;
        
        // 读取插值方法枚举值
        file >> n;
        if (n < INTERMAX)
            info->interMethod = (InterMethod)n;

        // 读取结束步数
        file >> n;
        info->endStep = n;

        // 读取输出时间间隔
        file >> nf;
        info->outputDt = nf;

        // 读取CFL数
        file >> nf;
        info->CFL = nf;

        // 读取案例编号
        file >> n;
        info->nCase = n;

        // 读取计算区域范围(x_min, x_max, y_min, y_max, z_min, z_max)
        real nf1, nf2, nf3, nf4, nf5, nf6;
        file >> nf1;
        file >> nf2;
        file >> nf3;
        file >> nf4;
        file >> nf5;
        file >> nf6;
        info->calZone = { nf1, nf2, nf3, nf4, nf5, nf6 };

        // 读取网格最大索引(i_max, j_max, k_max)
        int n1, n2, n3;
        file >> n1;
        file >> n2;
        file >> n3;
        info->iMax = { n1, n2, n3 };

        // 读取空间维度(1D/2D/3D)
        file >> n;
        info->dim = n;

        // 读取并设置OpenMP线程数
        file >> n;
        omp_set_num_threads(n);

        std::cout << "file mode initialization finished\n";

    } else {
        std::cout << "file mode initialization failed\n";
    }

    InterMethod interscheme;

    // ===================================求解器初始化与执行===================================

    // 初始化BlockSolver求解器对象
    BlockSolver bSolver(info);
    
    // 开始计时，记录求解器计算时间
    auto start = std::chrono::high_resolution_clock::now();
    
    // 根据方程类型选择时间推进方法
    // EULER方程使用基于CFL条件的自适应时间步长推进方法(stepsLoopCFL)
    // 其他方程使用固定时间步长推进方法(stepsLoop)
    if (info->eqType != EULER){
        bSolver.stepsLoop();
        // printf("if的");
        // std::cout << "if的" << '\n';
    }
    else{
        bSolver.stepsLoopCFL();
        // printf("else的");
    }
        // bSolver.Test();
        
    // 结束计时并计算总运行时间
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
                        .count();

    // bSolver.stepsLoopDTS();
    // bSolver.solve();
    // bSolver.outputPrim();
    // bSolver.Test();

    // ===================================结果输出与统计===================================

    // 创建timeInfo.txt文件用于输出计算时间统计信息
    std::ofstream timeinfo("timeInfo.txt");

    // 输出计算时间和步数信息到控制台
    std::cout << "totaltime= " << duration << "   Finish\n";
    std::cout << "time= " << timepp / 1e6 << "   Finish\n";
    std::cout << "timesteps= " << bSolver.timesteps << "   Finish\n";
    std::cout << "solvertime= " << timesss << '\n';

    // 输出计算时间和步数信息到timeInfo.txt文件
    timeinfo << info->interMethod << std::endl;
    timeinfo << "totaltime= " << duration << "   Finish\n";
    timeinfo << "time= " << timepp / 1e6 << "   Finish\n";
    timeinfo << "timesteps= " << bSolver.timesteps << "   Finish\n";
    timeinfo << "solvertime= " << timesss << '\n';
}