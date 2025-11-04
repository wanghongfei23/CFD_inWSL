
#include "blockSolver.hpp"
#include "eigenSystem.hpp"
#include <fstream>

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
    // info->interMethod = TCNS5;
    // info->interMethod=LINEAR5;
    // info->interMethod=WCNSZ5Char;
    // info->BVD=true;
    // info->interMethod = NICEST5;
    // info->interMethod=WCNS5CONG;
    //  info->interMethod=TCNSCongA;
    // info->interMethod = WCNS5CONGZ;
    // info->sourceType=GRAVITY;

    // 【王鸿飞】begin-1
    // info->interMethod = WHFTCNSA;
    info->interMethod = WHFTCNSAF002;
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
    // info->iMax = { 401, 401, 2 };
    // info->dim = 2;

    // Riemann 2 vortex
    //  info->endStep=1;
    //  info->outputDt=0.3;
    //  info->CFL=0.5;
    //  info->nCase=1;
    //  info->calZone={-0.5,0.5,-0.5,0.5,0,0};
    //  info->iMax={801,801,2};
    //  info->dim=2;

    // RT instability
    // 记得改GAMMA
     info->endStep=1;
     info->outputDt=1.95;
     info->CFL=0.5;
     info->nCase=3;
     info->calZone={0,0.25,0,1,0,0};
    //  info->iMax={201,801,2};
     info->iMax={65,257,2};
     info->dim=2;
     info->sourceType=GRAVITY;

    // info->diffMethod=HDS6;
    // Double Mach
    //  info->endStep=20;
    //  info->outputDt=0.01;
    //  info->CFL=0.5;
    //  info->nCase=4;
    //  info->calZone={0,4,0,1,0,0};
    //  info->iMax={801,201,2};
    //  info->dim=2;

    // file config mode
    std::ifstream file("info.txt");
    if (file.is_open()) {
        int n;
        real nf;
        file >> n;
        if (n < INTERMAX)
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

    std::ofstream timeinfo("timeInfo.txt");

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