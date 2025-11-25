
#include "blockSolver.hpp"
#include "eigenSystem.hpp"
#include <fstream>
// 【王鸿飞】begin-1命名
#include <map>
#include <string>
// 定义静态映射表，将算例映射到字符串
static std::map<int,std::string> exampleStr1D={
    {0,"Sod"},
    {1,"ShuOsher"},
    {2,"Lax"},
    {3,"sedov"},
    {4,"Woodward_Colella"},
    {5,"Double_sparse_wave"}
};

static std::map<int,std::string> exampleStr2D={
    {0,"2D_Riemann_1"},
    {1,"2D_Riemann_2"},
    {2,"implosion"},
    {3,"RTI"},
    {4,"Double_Mach"},
    {5,"2D_Riemann_3"},
    {6,"KHI"}
};

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
    {WHFTCNSASF102_par_01,"TENO-AS-myF102_par_01"},
    {WHFTCNSASF103_par_01,"TENO-AS-myF103_par_01"},
    {WHFTCNSAS_fx,"TENO-AS-fx"},
    {WHFTCNSAS_initial,"TENO-AS-initial"},
    {WHFTCNSAS_approx,"TENO-AS-approx"},
    {WHFTCNSAS_fx_real,"TENO-AS-fx-real"},
    {WHFTCNSAS_floor,"TENO-AS-floor"},
    {temp008,"temp_name_008"},
    {temp009,"temp_name_009"},
    {temp010,"temp_name_010"}
};
// 【王鸿飞】end-1命名

int main()
{
    // std::array<real,4> prim0={1.0,0.75,-0.5,1.0};
    // std::array<real,4> prim={2.0,-0.75,0.5,1.0};
    // eigensystemEuler2D eig=eigensystemEuler2D(prim0,{1,0,0});
    // auto eigValues=eig.primToChar(prim);
    // auto prim2=eig.charToPrim(eigValues);
    // std::cout<<"finish\n";

    // omp_set_num_threads(10);
    omp_set_num_threads(1);

    Info* info = new Info;

    info->eqType = EULER;
    info->spMethod = WCNS5;
    info->diffMethod = MND6;


    // info->diffMethod = TRAD6;
    // info->interMethod=LINEAR5;
    // info->interMethod=WCNSZ5Char;
    // info->BVD=true;
    // info->interMethod = NICEST5;
    // info->interMethod=WCNS5CONG;
    // info->sourceType=GRAVITY;

    // 【王鸿飞】begin-1

    // F102中的1表示改进版本 

    // info->interMethod = WCNS5; //weno5_JSchen
    // info->interMethod = WCNSZ5; //weno5_Z
    // info->interMethod = TCNS5; //Teno5_Z
    // info->interMethod = WCNS5CONGZ;//Teno5_CongZ
    // info->interMethod = WHFTCNSA;
    // info->interMethod = WHFTCNSASF002;
    // info->interMethod = WHFTCNSAH002;
    // info->interMethod = WHFTCNSASF102;
    // info->interMethod = WHFTCNSASF103;
    // info->interMethod = WHFTCNSASF102_par_01;
    // info->interMethod = WHFTCNSASF103_par_01;
    // info->interMethod = WHFTCNSAS_fx; // 不进行函数拟合，直接用原来近似的指数形式CT'
    // info->interMethod = WHFTCNSAS_initial; // 原封不动的叠加A和S
    info->interMethod = WHFTCNSAS_approx;
    // info->interMethod = WHFTCNSAS_fx_real;
    // info->interMethod = WHFTCNSAS_floor;
    // info->interMethod = temp008;
    // info->interMethod = temp009;
    // info->interMethod = temp010;

    // 【王鸿飞】end

    // Shu-Osher
    //  info->endStep=1;
    //  info->CFL=0.5;
    //  info->outputDt=1.8;
    //  info->nCase=1;
    //  info->calZone={0,10.0,0,0,0,0};
    //  info->iMax={201,2,2};
    //  info->dim=1;

    // // sod tube
    // info->CFL = 0.5;
    // info->endStep = 20;
    // info->outputDt = 0.01;
    // info->nCase = 0;
    // info->calZone = { -0.5, 0.5, 0, 0, 0, 0 };
    // info->iMax = { 201, 2, 2 };
    // info->dim = 1;

    // lax sod tube
    // info->endStep = 14;
    // info->outputDt = 0.01;
    // info->CFL = 0.5;
    // info->nCase = 2;
    // info->calZone = { -0.5, 0.5, 0, 0, 0, 0 };
    // info->iMax = { 201, 2, 2 };
    // info->dim = 1;

    // lax sod tube speed test
    //  info->endStep=14;
    //  info->outputDt=0.01;
    //  info->CFL=0.1;
    //  info->nCase=2;
    //  info->calZone={-0.5,0.5,0,0,0,0};
    //  info->iMax={2001,2,2};
    //  info->dim=1;

    // sedov
    //  info->endStep=1;
    //  info->outputDt=0.001;
    //  info->CFL=0.5;
    //  info->nCase=3;
    //  info->calZone={-2,2,0,0,0,0};
    //  info->iMax={400,2,2};
    //  info->dim=1;

    // Woodward-Colella
    // info->endStep = 1;
    // info->outputDt = 0.038;
    // info->CFL = 0.1;
    // info->nCase = 4;
    // info->calZone = { 0, 1, 0, 0, 0, 0 };
    // info->iMax = { 401, 2, 2 };
    // info->dim = 1;

    // 双稀疏波
    //  info->endStep=100;
    //  info->outputDt=0.01;
    //  info->CFL=0.5;
    //  info->nCase=5;
    //  info->calZone={-5,5,0,0,0,0};
    //  info->iMax={401,2,2};
    //  info->dim=1;

    // implosion
    //  info->endStep=25;
    //  info->outputDt=0.1;
    //  info->CFL=0.5;
    //  info->nCase=2;
    //  info->calZone={-0.3,0.3,-0.3,0.3,0,0};
    //  info->iMax={401,401,2};
    //  info->dim=2;

    // Riemann 1
    // info->endStep = 1;
    // info->outputDt = 0.8;
    // info->CFL = 0.5;
    // info->nCase = 0;
    // info->calZone = { -0.5, 0.5, -0.5, 0.5, 0, 0 };
    // info->iMax = { 401, 401, 2 };//参考
    // info->dim = 2;

    // Riemann 2 vortex
    //  info->endStep=1;
    //  info->outputDt=0.3;
    //  info->CFL=0.5;
    //  info->nCase=1;
    //  info->calZone={-0.5,0.5,-0.5,0.5,0,0};
    //  info->iMax={801,801,2};//参考
    //  info->dim=2;

    // Riemann 3
    //  info->endStep=8;
    //  info->outputDt=0.05;
    //  info->CFL=0.5;
    //  info->nCase=5;
    //  info->calZone={-0.5,0.5,-0.5,0.5,0,0};
    //  info->iMax={801,801,2};//参考
    //  info->dim=2;

    // RT instability
    // 记得改GAMMA
     info->endStep=1;
     info->outputDt=1.95;
     info->CFL=0.5;
     info->nCase=3;
     info->calZone={0,0.25,0,1,0,0};
    //  info->iMax={201,801,2};
     info->iMax={101,401,2};//参考
    //  info->iMax={65,257,2};
     info->dim=2;
     info->sourceType=GRAVITY;

    // info->diffMethod=HDS6;
    // Double Mach
    //  info->endStep=20;
    //  info->outputDt=0.01;
    //  info->CFL=0.5;
    //  info->nCase=4;
    //  info->calZone={0,4,0,1,0,0};
    //  info->iMax={801,201,2};//参考
    //  info->dim=2;

    // file config mode
    std::ifstream file("info.txt");
    if (file.is_open()) {
        int n;
        real nf;
        file >> n;
        if (n < temp010)
            info->interMethod = (InterMethod)n;

        file >> n;
        info->endStep = n;

        file >> nf;
        info->outputDt = nf;

        file >> nf;
        info->CFL = nf;

        file >> n;
        info->nCase = n;

        real nf1, nf2, nf3, nf4, nf5, nf6;
        file >> nf1;
        file >> nf2;
        file >> nf3;
        file >> nf4;
        file >> nf5;
        file >> nf6;
        info->calZone = { nf1, nf2, nf3, nf4, nf5, nf6 };

        int n1, n2, n3;
        file >> n1;
        file >> n2;
        file >> n3;
        info->iMax = { n1, n2, n3 };

        file >> n;
        info->dim = n;

        file >> n;
        omp_set_num_threads(n);

        std::cout << "file mode initialization finished\n";

    } else {
        std::cout << "file mode initialization failed\n";
    }

    InterMethod interscheme;

    BlockSolver bSolver(info);
    auto start = std::chrono::high_resolution_clock::now();
    if (info->eqType != EULER)
        bSolver.stepsLoop();
    else
        bSolver.stepsLoopCFL();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    // bSolver.stepsLoopDTS();
    // bSolver.solve();
    // bSolver.outputPrim();
    // bSolver.Test();

    // 原命名
    // std::ofstream timeinfo("timeInfo.txt");

    // 【王鸿飞】begin-1命名
    // 创建与info.cpp中filename()函数类似命名方式的文件名
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
    // 【王鸿飞】end-1命名
    
    std::cout << "totaltime= " << duration << "   Finish\n";
    std::cout << "time= " << timepp / 1e6 << "   Finish\n";
    std::cout << "timesteps= " << bSolver.timesteps << "   Finish\n";
    std::cout << "solvertime= " << timesss << '\n';

    timeinfo << info->interMethod << std::endl;
    timeinfo << "totaltime= " << duration << "   Finish\n";
    timeinfo << "time= " << timepp / 1e6 << "   Finish\n";
    timeinfo << "timesteps= " << bSolver.timesteps << "   Finish\n";
    timeinfo << "solvertime= " << timesss << '\n';
}