#include "differ.hpp"
#include "differenceScheme.hpp"
#include <boost/circular_buffer.hpp>
#include <ranges>
#include <span>

void Differ::solve()
{
    assert(rhsR.getN() + 1 == fluxesHalf->getN());
    int n = rhsR.getN(), nvar = rhsR.getNVar();
    real h = constantH;
    auto iterFluxHalf = fluxesHalf->begin();
    auto iterFluxNode = fluxesNode->begin();

    for (auto iterRhs = rhsR(0); iterRhs < rhsR(n); iterRhs += rhsR.getOffset()) {

        for (int ivar = 0; ivar < nvar; ivar++) {
            iterRhs[ivar] += secondOrder(iterFluxHalf + ivar, iterFluxNode + ivar, nvar) / h;
        }
    }
}

void ExplicitDif6::solve()
{
    assert(rhsR.getN() + 5 == fluxesHalf->getN());

    int n = rhsR.getN(), nvar = rhsR.getNVar();
    real h = constantH;
    auto iterFluxHalf = fluxesHalf->begin();
    constexpr std::array<real, 3> w = { 75.0 / 64.0, 25.0 / (128.0 * 3.0),
        3.0 / (128.0 * 5.0) };
    constexpr int nStencil = 6;
    boost::circular_buffer<std::span<real>> varArrays(nStencil);
    for (int i : std::views::iota(0, nStencil - 1)) {
        varArrays.push_back(std::span<real>(iterFluxHalf, nvar));
        iterFluxHalf += nvar;
    }
    assert(varArrays.size() == nStencil - 1);

    for (auto iterRhs = rhsR(0); iterRhs < rhsR(n); iterRhs += rhsR.getOffset()) {
        varArrays.push_back(std::span<real>(iterFluxHalf, nvar));
        for (int ivar = 0; ivar < nvar; ivar++) {
            iterRhs[ivar] += (w[0] * (varArrays[3][ivar] - varArrays[2][ivar])
                                 - w[1] * (varArrays[4][ivar] - varArrays[1][ivar])
                                 + w[2] * (varArrays[5][ivar] - varArrays[0][ivar]))
                / h;
        }

        iterFluxHalf += nvar;
    }
    assert(iterFluxHalf == fluxesHalf->end());
}

void MidNodeAndNodeDif6::solve()
{
    assert(rhsR.getN() + 3 == fluxesHalf->getN());
    assert(rhsR.getN() + 2 == fluxesNode->getN());
    int n = rhsR.getN(), nvar = rhsR.getNVar();
    real h = constantH;
    constexpr std::array<real, 3> w = { 3.0 / 2.0, -3.0 / 10.0, 1.0 / 30.0 };
    auto iterFluxHalf = fluxesHalf->begin();
    auto iterFluxNode = fluxesNode->begin();

    // 初始化数组
    constexpr int nStencilHalf = 4;
    constexpr int nStencilNode = 3;
    boost::circular_buffer<std::span<real>> varArraysHalf(nStencilHalf),
        varArraysNode(nStencilNode);
    for (int i : std::views::iota(0, nStencilHalf - 1)) {
        varArraysHalf.push_back(std::span<real>(iterFluxHalf, nvar));
        iterFluxHalf += nvar;
    }
    for (int i : std::views::iota(0, nStencilNode - 1)) {
        varArraysNode.push_back(std::span<real>(iterFluxNode, nvar));
        iterFluxNode += nvar;
    }
    assert(varArraysHalf.size() == nStencilHalf - 1);
    assert(varArraysNode.size() == nStencilNode - 1);

    auto iterRhs = rhsR(0);
    for (; iterRhs < rhsR(n); iterRhs += rhsR.getOffset()) {
        varArraysHalf.push_back(std::span<real>(iterFluxHalf, nvar));
        varArraysNode.push_back(std::span<real>(iterFluxNode, nvar));
        for (int ivar : std::views::iota(0, nvar)) {
            //chenyuqing:差分的实际逻辑
            iterRhs[ivar] += (w[1] * (varArraysNode[2][ivar] - varArraysNode[0][ivar])
                                 + w[0] * (varArraysHalf[2][ivar] - varArraysHalf[1][ivar])
                                 + w[2] * (varArraysHalf[3][ivar] - varArraysHalf[0][ivar]))
                / h;
            //chenyuqing：修改差分公式
                // iterRhs[ivar] += w[1]*varArraysNode[2][ivar] / h - w[1]*varArraysNode[0][ivar] / h
                //                  - w[0]* varArraysHalf[1][ivar] / h + w[0] *varArraysHalf[2][ivar] / h 
                //                  + w[2] * varArraysHalf[3][ivar] / h - w[2] * varArraysHalf[0][ivar]/ h
                // ;

        }

        iterFluxHalf += nvar;
        iterFluxNode += nvar;
    }
    assert(iterRhs == rhsR(n));
    assert(iterFluxHalf == fluxesHalf->end());
}