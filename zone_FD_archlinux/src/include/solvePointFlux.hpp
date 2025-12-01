#pragma once
#include "dataManipulater.hpp"
#include <concepts>

/*------------------concepts bigin---------------------------*/
template <typename T, std::size_t NVar>
concept Flux = requires(T f, std::span<real, NVar> arr, std::span<real, NVar> arr2,
    std::array<real, 3> norm) { T(arr, arr2, norm); };

template <std::size_t NVar>
using FluxN = void (*)(const std::span<real, NVar>&,
    const std::span<real, NVar>&,
    const std::array<real, 3>&);
/*-----------------concepts end--------------------------------*/

template <std::size_t NVar, FluxN<NVar> fluxSolver>
class SolvePointSolver : public DataManipulater {
public:
    SolvePointSolver() {};
    SolvePointSolver(DataReader varsReader_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_, int nl_, int nr_)
        : DataManipulater(varsReader_, bndL_, bndR_)
        , nl(nl_)
        , nr(nr_) {};

    void setNlNr(int nl_, int nr_)
    {
        nl = nl_;
        nr = nr_;
    }
    void check() { std::cout << "initialized successfully SolvePointSolver\n"; }

    void solve() final
    {
        if (nl == 0 && nr == 0)
            return;
        initVarsR();
        iter = data->begin();
        leftBnd();
        internal();
        rightBnd();
        // std::cout << std::endl;
        assert(iter == data->end());
    };
    void initVarsR() final;
    void leftBnd() final;
    void internal() final;
    void rightBnd() final;
    void solveI(std::vector<real>::iterator);

protected:
    int nl = 0, nr = 0;
};
template <std::size_t NVar, FluxN<NVar> fluxSolver>
void SolvePointSolver<NVar, fluxSolver>::initVarsR()
{
    if (!data) {
        data = std::make_shared<Data>(n + nl + nr, NVar);
    }
    assert(NVar == nvar);
}

template <std::size_t NVar, FluxN<NVar> fluxSolver>
void SolvePointSolver<NVar, fluxSolver>::leftBnd()
{
    for (int i = nl - 1; i >= 0; i--) {
        auto iterBnd = (*bndL)(i);
        solveI(iterBnd);
    }
}

template <std::size_t NVar, FluxN<NVar> fluxSolver>
void SolvePointSolver<NVar, fluxSolver>::internal()
{
    auto endIter = varsReader(n);
    for (auto iterVar = varsReader(0); iterVar != endIter; iterVar += varsReader.getOffset()) {
        solveI(iterVar);
        assert(iterVar <= endIter);
    }
}

template <std::size_t NVar, FluxN<NVar> fluxSolver>
void SolvePointSolver<NVar, fluxSolver>::rightBnd()
{
    for (int i = 0; i < nr; i++) {
        auto iterBnd = (*bndR)(i);
        solveI(iterBnd);
    }
}

template <std::size_t NVar, FluxN<NVar> fluxSolver>
void SolvePointSolver<NVar, fluxSolver>::solveI(
    std::vector<real>::iterator iterBnd)
{
    assert(NVar == nvar);
    std::span<real, NVar> input(iterBnd, NVar), output(iter, NVar);
    (*fluxSolver)(input, output, norm);
    iter += NVar;
    // std::cout << iter - data->end() << '\n';
}