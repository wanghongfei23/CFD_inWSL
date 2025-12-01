#pragma once
#include "eigenSystem.hpp"
#include "interScheme.hpp"
#include "reconstructor.hpp"
#include <functional>
#include "positive_preserving.hpp"

typedef real (*InterpolationScheme5Order)(std::array<real, 5>);

template <InterpolationScheme5Order inter5>
class Recon5OrderFaceCenter : public Reconstuctor {
public:
    Recon5OrderFaceCenter() {};

    Recon5OrderFaceCenter(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Reconstuctor(varsR_, bndL_, bndR_) {};

protected:
    void initVarsR() override;
    void leftBnd() override;
    void internal() override;
    void rightBnd() override;
    virtual void reconI(std::array<std::vector<real>::iterator, 6>);
    int tt = 0;
};

template <InterpolationScheme5Order inter5>
class Recon5Order1DEulerEig : public Recon5OrderFaceCenter<inter5> {

public:
    Recon5Order1DEulerEig() {};
    Recon5Order1DEulerEig(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Recon5OrderFaceCenter<inter5>(varsR_, bndL_, bndR_) {};

protected:
    const int NVAR = 3;
    void reconI(std::array<std::vector<real>::iterator, 6>) final;
};

template <InterpolationScheme5Order inter5>
class Recon5Order2DEulerEig : public Recon5OrderFaceCenter<inter5> {
public:
    Recon5Order2DEulerEig() {};
    Recon5Order2DEulerEig(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : Recon5OrderFaceCenter<inter5>(varsR_, bndL_, bndR_) {};

protected:
    const int NVAR = 4;
    void reconI(std::array<std::vector<real>::iterator, 6>) final;
};

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::initVarsR()
{
    if (!data) {
        int nPointR = bndL->getN() + n + bndR->getN() - 5;
        data = std::make_shared<Data>(nPointR, nvar * 2);
    }
}

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::leftBnd()
{
    // std::cout << (iter - data->end()) / 2 / nvar << '\n';
    int nBnd = bndL->getN();
    std::vector<std::vector<real>::iterator> tempIters(nBnd + 5);
    auto tempItersIter = tempIters.begin();
    for (int i = 0; i < nBnd; i++) {
        int iInver = nBnd - 1 - i;
        auto tempppp = (*bndL)(iInver);
        (*tempItersIter++) = (*bndL)(iInver);
    }
    for (int i = 0; i < 5; i++) {
        (*tempItersIter++) = varsReader(i);
    }
    assert(tempItersIter == tempIters.end());

    for (int i = 0; i < nBnd; i++) {
        reconI({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4], tempIters[i + 5] });
    }
}

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::internal()
{

    for (int i = 0; i < n - 5; i++) {
        auto iter = varsReader(i);
        reconI({
            varsReader(i + 0),
            varsReader(i + 1),
            varsReader(i + 2),
            varsReader(i + 3),
            varsReader(i + 4),
            varsReader(i + 5),
        });
    }
}

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::rightBnd()
{
    // std::cout << (iter - data->end()) / nvar / 2 << '\n';
    int nBnd = bndR->getN();
    std::vector<std::vector<real>::iterator> tempIters(nBnd + 5);
    auto tempItersIter = tempIters.begin();
    for (int i = 0; i < 5; i++) {
        (*tempItersIter++) = varsReader(n - 5 + i);
    }
    for (int i = 0; i < nBnd; i++) {
        (*tempItersIter++) = (*bndR)(i);
    }

    assert(tempItersIter == tempIters.end());
    for (int i = 0; i < nBnd; i++) {
        reconI({ tempIters[i + 0], tempIters[i + 1], tempIters[i + 2],
            tempIters[i + 3], tempIters[i + 4], tempIters[i + 5] });
    }

    assert(iter == data->end());
}

template <InterpolationScheme5Order inter5>
void Recon5OrderFaceCenter<inter5>::reconI(
    std::array<std::vector<real>::iterator, 6> input)
{
    // 每一次运行iter都应该加2*nVar哟
    // 这里实现了原始变量重构
    for (int i = 0; i < nvar; i++) {
        *iter = inter5({
            input[0][i],
            input[1][i],
            input[2][i],
            input[3][i],
            input[4][i],
        });
        *(iter + nvar) = inter5({
            input[5][i],
            input[4][i],
            input[3][i],
            input[2][i],
            input[1][i],
        });
        iter++;
    }
    iter += nvar;
}

template <InterpolationScheme5Order inter5>
void Recon5Order1DEulerEig<inter5>::reconI(
    std::array<std::vector<real>::iterator, 6> input)
{
    // 每一次运行iter都应该加2*nVar哟
    enum { R,
        U,
        P };
    auto primL = input[2], primR = input[3];
    eigensystemEuler1D eig = eigensystemEuler1D({ primL[R], primL[U], primL[P] },
        { primR[R], primR[U], primR[P] });
    std::array<real, 5> q1L, q2L, q3L, q1R, q2R, q3R;
    for (int j = 0; j < 6; j++) {

        auto charTemp = eig.primToChar({ input[j][R], input[j][U], input[j][P] });

        if (j < 5) {
            q1L[j] = charTemp[0];
            q2L[j] = charTemp[1];
            q3L[j] = charTemp[2];
        }

        int jj = 5 - j;
        if (jj < 5) {
            q1R[jj] = charTemp[0];
            q2R[jj] = charTemp[1];
            q3R[jj] = charTemp[2];
        }
    }
    // auto start = std::chrono::steady_clock::now();
    auto Q1LL = inter5(q1L);
    auto Q1RR = inter5(q1R);
    auto Q2LL = inter5(q2L);
    auto Q2RR = inter5(q2R);
    auto Q3LL = inter5(q3L);
    auto Q3RR = inter5(q3R);
    // auto stop = std::chrono::steady_clock::now();
    // auto duration =
    // std::chrono::duration_cast<std::chrono::nanoseconds>(stop
    // - start).count(); timep+=duration;

    // auto resTempL = eig.charToPrim({ Q1LL, Q2LL, Q3LL });
    // auto resTempR = eig.charToPrim({ Q1RR, Q2RR, Q3RR });
    // if (resTempL[0]<0||resTempL[2]<0) resTempL={ primL[R], primL[U], primL[P] };
    // if (resTempR[0]<0||resTempR[2]<0) resTempR={ primR[R], primR[U], primR[P] };

// 将原始变量（密度、速度、压力）转换为守恒变量
auto consL = eig.primToCons({ primL[R], primL[U], primL[P] });
auto consR = eig.primToCons({ primR[R], primR[U], primR[P] });

// 定义 lambda 表达式，提取守恒变量中的密度
auto funcRho = [](std::array<real, 3> const& W) { return W[0]; };

// 定义 lambda 表达式，计算密度乘以声速的平方
auto funcRhoCsq = [](std::array<real, 3> const& W) {
    real rho = W[0];          // 密度
    real rhou = W[1];         // 动量
    real rhoE = W[2];         // 总能

    // 计算速度
    real u = rhou / rho;      // 速度
    // 计算动能
    real ekin = 0.5 * u * u;  // 动能
    // 计算压力
    real p = (GAMMA - 1) * (rhoE - rho * ekin); // 压力
    // 返回密度乘以声速的平方
    return GAMMA * p;
};

// 将特征变量转换为守恒变量
auto consTempL = eig.charToCons({ Q1LL, Q2LL, Q3LL });
auto consTempR = eig.charToCons({ Q1RR, Q2RR, Q3RR });

// 对守恒变量进行正限制，确保密度为正
variable_positive_limiter<3>(consL, consTempL, funcRho);
variable_positive_limiter<3>(consR, consTempR, funcRho);

// 对守恒变量进行正限制，确保密度乘以声速的平方为正
variable_positive_limiter<3>(consL, consTempL, funcRhoCsq);
variable_positive_limiter<3>(consR, consTempR, funcRhoCsq);

// 将守恒变量转换回原始变量
auto resTempL = eig.consToPrim(consTempL);
auto resTempR = eig.consToPrim(consTempR);

// 将结果复制到指定位置
std::copy(resTempL.begin(), resTempL.end(), this->iter);
std::copy(resTempR.begin(), resTempR.end(), this->iter + this->nvar);

// 更新迭代器
this->iter += 2 * this->nvar;
}

template <InterpolationScheme5Order inter5>
void Recon5Order2DEulerEig<inter5>::reconI(
    std::array<std::vector<real>::iterator, 6> input)
{
    // 每一次运行iter都应该加2*nVar哟
    //chenyuqing: 重构的业务逻辑
    assert(this->nvar == 4);
    enum { R,
        U,
        V,
        P };
    auto primL = input[2], primR = input[3];
    //chenyuqing: 创建特征变换对象eigensystemEuler2D
    eigensystemEuler2D eig = eigensystemEuler2D({ primL[R], primL[U], primL[V], primL[P] },
        { primR[R], primR[U], primR[V], primR[P] }, this->norm);
    std::array<real, 5> q1L, q2L, q3L, q4L, q1R, q2R, q3R, q4R;
    for (int j = 0; j < 6; j++) {
        //chenyuqing: 把原始变量转化为特征变量
        auto charTemp = eig.primToChar({ input[j][R], input[j][U], input[j][V], input[j][P] });

        //chenyuqing: 把特征变量转化为合适的数组，包含了模板翻转
        if (j < 5) {
            q1L[j] = charTemp[R];
            q2L[j] = charTemp[U];
            q3L[j] = charTemp[V];
            q4L[j] = charTemp[P];
        }

        int jj = 5 - j;
        if (jj < 5) {
            q1R[jj] = charTemp[R];
            q2R[jj] = charTemp[U];
            q3R[jj] = charTemp[V];
            q4R[jj] = charTemp[P];
        }
    }
    // auto start = std::chrono::steady_clock::now();
    //chenyuqing: 进行非线性插值
    auto Q1LL = inter5(q1L);
    auto Q1RR = inter5(q1R);
    auto Q2LL = inter5(q2L);
    auto Q2RR = inter5(q2R);
    auto Q3LL = inter5(q3L);
    auto Q3RR = inter5(q3R);
    auto Q4LL = inter5(q4L);
    auto Q4RR = inter5(q4R);
    // auto stop = std::chrono::steady_clock::now();
    // auto duration =
    // std::chrono::duration_cast<std::chrono::nanoseconds>(stop
    // - start).count(); timep+=duration;
    //chenyuqing: 变回原始变量
    auto resTempL = eig.charToPrim({ Q1LL, Q2LL, Q3LL, Q4LL });
    auto resTempR = eig.charToPrim({ Q1RR, Q2RR, Q3RR, Q4RR });

    // //chenyuqing: 简单的保正限制器
    if (resTempL[0]<0||resTempL[3]<0) resTempL={ primL[R], primL[U], primL[V], primL[P] };
    if (resTempR[0]<0||resTempR[3]<0) resTempR={ primR[R], primL[U], primL[V], primR[P] };

    // 将原始变量（密度、速度分量 u 和 v、压力）转换为守恒变量
// auto consL = eig.primToCons({ primL[R], primL[U], primL[V], primL[P] });
// auto consR = eig.primToCons({ primR[R], primR[U], primR[V], primR[P] });

// 定义 lambda 表达式，提取守恒变量中的密度
// auto funcRho = [](std::array<real, 4> const& W) { return W[0]; };

// // 定义 lambda 表达式，计算密度乘以声速的平方
// auto funcRhoCsq = [](std::array<real, 4> const& W) {
//     real rho = W[0];          // 密度
//     real rhou = W[1];         // x 方向动量
//     real rhov = W[2];         // y 方向动量
//     real rhoE = W[3];         // 总能

//     // 计算速度分量
//     real u = rhou / rho;      // x 方向速度
//     real v = rhov / rho;      // y 方向速度
//     // 计算动能
//     real ekin = 0.5 * (u * u + v * v);  // 动能
//     // 计算压力
//     real p = (GAMMA - 1) * (rhoE - rho * ekin); // 压力
//     // 返回密度乘以声速的平方
//     return GAMMA * p;
// };

// // 将特征变量转换为守恒变量
// auto consTempL = eig.charToCons({ Q1LL, Q2LL, Q3LL, Q4LL });
// auto consTempR = eig.charToCons({ Q1RR, Q2RR, Q3RR, Q4RR });

// // 对守恒变量进行正限制，确保密度为正
// variable_positive_limiter<4>(consL, consTempL, funcRho);
// variable_positive_limiter<4>(consR, consTempR, funcRho);

// // 对守恒变量进行正限制，确保密度乘以声速的平方为正
// variable_positive_limiter<4>(consL, consTempL, funcRhoCsq);
// variable_positive_limiter<4>(consR, consTempR, funcRhoCsq);

// // 将守恒变量转换回原始变量
// auto resTempL = eig.consToPrim(consTempL);
// auto resTempR = eig.consToPrim(consTempR);

// 将结果复制到指定位置
std::copy(resTempL.begin(), resTempL.end(), this->iter);
std::copy(resTempR.begin(), resTempR.end(), this->iter + this->nvar);

// 更新迭代器
this->iter += 2 * this->nvar;

    // //chenyuqing: 把重构后的左右值放到总的数组里
    // std::copy(resTempL.begin(), resTempL.end(), this->iter);
    // std::copy(resTempR.begin(), resTempR.end(), this->iter + this->nvar);
    // this->iter += 2 * this->nvar;
}