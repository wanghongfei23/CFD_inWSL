#pragma once

#include <concepts>
#include <span>
#include <type_traits>

#include "data.hpp"
#include "info.hpp"

/*------------------concepts bigin---------------------------*/

template <std::size_t NVar>
using RiemannSolverN = void (*)(const std::span<real, NVar * 2>&,
    const std::span<real, NVar>&,
    const std::array<real, 3>&);

template <typename T, std::size_t NVar>
concept RiemannSolver = requires(T f, std::span<real, NVar> arr, std::array<real, 3> norm) {
    { T(arr, norm) } -> std::same_as<std::array<real, NVar>>;
};

/*-----------------concepts end--------------------------------*/

// 通量点求解器模板类，用于计算通量点上的数值通量
template <std::size_t NVar, RiemannSolverN<NVar> solver>
class FluxPointSolver {
public:
    // 构造函数
    FluxPointSolver() {};
    FluxPointSolver(std::shared_ptr<Data> valsR_, int idim_, int nvar_)
        : valsR(valsR_)
        , idim(idim_)
        , nvar(nvar_) {};

    // 初始化函数
    void init(std::shared_ptr<Data> valsR_, int idim_, int nvar_)
    {
        valsR = valsR_;
        idim = idim_;
        nvar = nvar_;
    }
    
    // 数据获取函数
    std::shared_ptr<Data> getData() { return fluxes; }
    
    // 设置法向量函数
    void setConstNorm(std::array<real, 3> norm_);
    
    // 检查函数
    void check() { std::cout << "initialized successfully FluxPointSolver\n"; }
    
    // 求解函数
    void solve();
    
    // 通量数据指针
    std::shared_ptr<Data> fluxes;

protected:
    // 数据和参数
    std::shared_ptr<Data> valsR;      // 值数据指针
    int idim, nvar;                   // 维度索引和变量数
    std::array<real, 3> norm;         // 法向量
};

// 设置常数法向量实现
template <std::size_t NVar, RiemannSolverN<NVar> solver>
void FluxPointSolver<NVar, solver>::setConstNorm(std::array<real, 3> norm_)
{
    norm = norm_;
}

// 求解实现
template <std::size_t NVar, RiemannSolverN<NVar> solver>
void FluxPointSolver<NVar, solver>::solve()
{
    auto end = valsR->end();
    assert(nvar == NVar);
    assert(valsR->size() % 2 == 0);
    if (!fluxes)
        fluxes = std::make_shared<Data>(valsR->getN(), NVar);
    auto ivar = valsR->begin();
    auto iflux = fluxes->begin();
    for (; ivar != end; ivar += 2 * NVar, iflux += NVar) {
        std::span<real, 2 * NVar> input(ivar, 2 * NVar);
        std::span<real, NVar> output(iflux, NVar);
        solver(input, output, norm);
    }
    assert(ivar == end);
    assert(iflux == fluxes->end());
}