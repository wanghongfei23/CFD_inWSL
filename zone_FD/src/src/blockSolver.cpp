#include "blockSolver.hpp"
#include <fstream>

BlockSolver::BlockSolver()
{

    info = new Info();
    block = new Block();
    initer = new Initializer(info);
    eqn = new Equation();
    bnds = new Bnds();
    spDis = new SpDistributor();

    initer->initUniformBlock(block);
    initer->initEqution(eqn, block);
    initer->initBnds(bnds, eqn, block->getICMax(), block);
    initer->initSpDistributor(spDis, eqn, block, bnds);
    sourceTerm = new SourceTerm(eqn->getPrim(), eqn->getRhs(), info);

    cons = eqn->getCons();
    rhs = eqn->getRhs();
}

BlockSolver::BlockSolver(Info* info_)
{
    info = info_;
    block = new Block();
    initer = new Initializer(info);
    eqn = new Equation();
    bnds = new Bnds();
    spDis = new SpDistributor();

    initer->initUniformBlock(block);
    initer->initEqution(eqn, block);
    initer->initBnds(bnds, eqn, block->getICMax(), block);
    initer->initSpDistributor(spDis, eqn, block, bnds);
    sourceTerm = new SourceTerm(eqn->getPrim(), eqn->getRhs(), info);

    cons = eqn->getCons();
    rhs = eqn->getRhs();
}

void BlockSolver::RK3_SSP(real dt)
{

    Data tempdata(*cons);
    int n = cons->size();

    // third order RK
    // stage 1
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    sourceTerm->calSource();
// cgnsIO.BlockCgnsOutput(block,info);
// cgnsIO.solCgnsOutput(rhs,info);
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        (*cons)[i] = tempdata[i] - dt * (*rhs)[i];
    }
    info->t += dt;

    // stage 2
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    sourceTerm->calSource();
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        (*cons)[i] = 0.75 * tempdata[i] - 0.25 * dt * (*rhs)[i] + 0.25 * (*cons)[i];
    }
    info->t -= dt / 2;

    // stage 3
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    sourceTerm->calSource();
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        (*cons)[i] = 1.0 / 3.0 * tempdata[i] - 2.0 / 3.0 * dt * (*rhs)[i] + 2.0 / 3.0 * (*cons)[i];
    }
    info->t -= dt / 2;
}

void BlockSolver::DTS_Euler(real dt)
{

    int n = cons->size();
    Data tempdata(*cons);
    Data tempRhs(*cons);
    real maxres = 0;

    int imStep = 0;
    do {

        auto dtau = calLocalCFL();
        rhs->setZeros();
        eqn->consToPrim();
        bnds->update();
        spDis->rhsSolve();
        sourceTerm->calSource();
        // cgnsIO.BlockCgnsOutput(block,info);
        // cgnsIO.solCgnsOutput(rhs,info);
        tempRhs.setZeros();
#pragma omp parallel for
        for (int i = 0; i < n; i++) {
            real coef = dtau[i] / (dt + dtau[i]);
            tempRhs[i] = (tempdata[i] - (*cons)[i] - dt * (*rhs)[i]) * coef;
            (*cons)[i] = (*cons)[i] + tempRhs[i];
        }
        maxres = tempRhs.getLinf(0);
        std::cout << std::format("time = {:.4f} dt={:.10f} imStep={} rhoL2 = {}  \n", info->t, dt, imStep, maxres);
    } while (imStep++ < info->maxImplicitStep && maxres > 1e-7);
}

void BlockSolver::solve(real dt)
{
    RK3_SSP(dt);
    info->t += dt;
}

BlockSolver::~BlockSolver()
{
    delete info;
    delete block;
    delete initer;
    delete eqn;
    delete bnds;
    delete spDis;
    delete sourceTerm;
}

void BlockSolver::outputGrid()
{
    cgnsIO.BlockCgnsOutput(block, info);
}

void BlockSolver::outputCons()
{
    cgnsIO.solCgnsOutput(eqn->getCons(), info);
}
void BlockSolver::outputPrim()
{
    eqn->consToPrim();
    cgnsIO.solCgnsOutput(eqn->getPrim(), info);
}

void BlockSolver::stepsLoop()
{
    for (; info->step < info->endStep; info->step++) {
        if (info->step % info->outputInterval == 0) {
            outputGrid();
            outputPrim();
        }

        real dt = info->dt;
        solve(dt);
        timesteps++;
        std::cout << std::format("time = {:.4f}  rhoL2 = {:.4f}  \n", info->t, rhs->getL2(0));
    }
    outputGrid();
    outputPrim();
}
 
// 时间推进循环
void BlockSolver::stepsLoopCFL()
{
    // 定义一个极小值，用于浮点数比较
    const real eps = 1e-15;
    
    // 时间步进循环，直到达到结束步数
    while (info->step < info->endStep) {
        // 当达到输出时间时，输出网格和物理量信息
        if (info->outputT < eps) {
            outputGrid();
            outputPrim();
        }

        // 根据CFL条件计算时间步长
        auto dt = getTimeIntervalExplicit();
        // dt=info->dt;
        
        // 检查是否需要调整时间步长以匹配输出间隔
        if (info->outputT + dt > info->outputDt) {
            dt = info->outputDt - info->outputT;
            info->outputT = 0.0;
            info->step++;
        } else
            info->outputT += dt;

        // 执行一次时间步进求解并计时
        auto start = std::chrono::high_resolution_clock::now();
        solve(dt);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        timesss += duration;
        timesteps++;
        
        // 输出当前时间、时间步长和密度残差的L∞范数
        std::string output_str = std::format("time = {:.4f} dt={:.10f}  rhoLinf = {:.4f}  \n", info->t, dt, rhs->getLinf(0));
        
        std::cout << output_str;
        
        // 【王鸿飞】begin辅助工具
        
        // 将相同信息写入txt文件，使用与CGNS文件相同的命名方式，但不包含时间戳
        // std::string base_filename = info->filename();
        // // 移除时间戳部分（从最后一个" - t="开始的部分）
        // size_t pos = base_filename.rfind(" - t=");
        // if (pos != std::string::npos) {
        //     base_filename = base_filename.substr(0, pos) + ".txt";
        // } else {
        //     // 如果没有找到时间戳模式，则直接替换扩展名
        //     pos = base_filename.find(".cgns");
        //     if (pos != std::string::npos) {
        //         base_filename.replace(pos, 5, ".txt");
        //     }
        // }
        
        // std::ofstream log_file(base_filename, std::ios::app);
        // if (log_file.is_open()) {
        //     log_file << output_str;
        //     log_file.close();
        // }
        // 【王鸿飞】end辅助工具
    }
    
    // 循环结束后最后输出一次结果
    outputGrid();
    outputPrim();
}

void BlockSolver::Test()
{
    outputGrid();
    RK3_SSP(info->dt);

    // outputCons();
    // outputPrim();
    cgnsIO.solCgnsOutput(rhs, info);
}

std::vector<real> BlockSolver::calLocalCFL()
{

    real dt, lambda = 0;

    Data* prim = eqn->getPrim();
    int n = prim->getN();
    std::vector<real> res(n);
    eqn->consToPrim();
#pragma omp parallel for reduction(max : lambda)
    for (int i = 0; i < n; i++) {
        real iLambda;
        real dh = block->getMinDh(i);
        real RT = std::sqrt(GAMMA * (*prim)(i, 2) / (*prim)(i, 0)) + std::abs((*prim)(i, 1));
        if (info->dim == 2)
            iLambda = (std::sqrt(GAMMA * (*prim)(i, 3) / (*prim)(i, 0)) + std::abs((*prim)(i, 1)) + std::abs((*prim)(i, 2))) / dh;
        else if (info->dim == 1)
            iLambda = (std::sqrt(GAMMA * (*prim)(i, 2) / (*prim)(i, 0)) + std::abs((*prim)(i, 1))) / dh;
        res[i] = info->implicitCFL / iLambda;
    }
    if (info->eqType != EULER) {
        std::cout << "Block Solver error: CFL Loop with not EULER solver\n";
    }
    return res;
}

void BlockSolver::stepsLoopDTS()
{
    const real eps = 1e-15;
    while (info->step < info->endStep) {
        if (info->outputT < eps) {
            outputGrid();
            outputPrim();
        }

        real dt, lambda = 0;
        Data* prim = eqn->getPrim();
        int n = prim->getN();
        eqn->consToPrim();
#pragma omp parallel for reduction(max : lambda)
        for (int i = 0; i < n; i++) {
            real iLambda;
            real dh = block->getMinDh(i);
            real RT = std::sqrt(GAMMA * (*prim)(i, 2) / (*prim)(i, 0)) + std::abs((*prim)(i, 1));
            if (info->dim == 2)
                iLambda = (std::sqrt(GAMMA * (*prim)(i, 3) / (*prim)(i, 0)) + std::abs((*prim)(i, 1)) + std::abs((*prim)(i, 2))) / dh;
            else if (info->dim == 1)
                iLambda = (std::sqrt(GAMMA * (*prim)(i, 2) / (*prim)(i, 0)) + std::abs((*prim)(i, 1))) / dh;
            lambda = std::max(lambda, iLambda);
        }
        if (info->eqType != EULER) {
            std::cout << "Block Solver error: CFL Loop with not EULER solver\n";
            lambda = 1;
        }
        dt = info->CFL / lambda;
        // dt=info->dt;
        if (info->outputT + dt > info->outputDt) {
            dt = info->outputDt - info->outputT;
            info->outputT = 0.0;
            info->step++;
        } else
            info->outputT += dt;

        DTS_Euler(dt);
        info->t += dt;
        // std::cout<<std::format("time = {:.4f} dt={:.10f}  rhoL2 = {:.4f}  \n",info->t,dt,rhs->getL2(0));
    }
    outputGrid();
    outputPrim();
}

real BlockSolver::getTimeIntervalExplicit()
{
    real dt, lambda = 0;
    int dim = info->dim;
    if (info->eqType == EULER) {
        Data* prim = eqn->getPrim();
        int n = prim->getN();
        eqn->consToPrim();
        std::vector<real> lambdas(dim);
        for (int idim = 0; idim < dim; idim++) {
            real maxLambda = 0;
#pragma omp parallel for reduction(max : maxLambda)
            for (int i = 0; i < n; i++) {
                real iLambda;
                auto dhs = block->getCellInterval(i);
                real dh = info->constH ? info->geth(idim) : dhs.at(idim);
                iLambda = (std::sqrt(GAMMA * (*prim)(i, dim + 1) / (*prim)(i, 0)) + std::abs((*prim)(i, idim + 1))) / dh;
                maxLambda = std::max(maxLambda, iLambda);
            }
            lambdas[idim] = maxLambda;
        }
        lambda = 0;
        for (auto iLambda : lambdas)
            lambda += iLambda;

    } else {
        std::cout << "Block Solver error: CFL Loop with not EULER solver\n";
        lambda = 1;
    }
    dt = info->CFL / lambda;
    return dt;
}