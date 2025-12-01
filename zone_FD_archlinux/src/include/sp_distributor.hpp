#pragma once
#include "bnds.hpp"
#include "solverType.hpp"
#include <omp.h>

class SpDistributor {
public:
    void rhsSolve();

private:
    friend class Initializer;

    int nCons, nPrim, dim;
    std::array<int, 3> iMax;
    Data *prim, *cons;
    Data* rhs;

    Bnds* bnds;
    Info* info;
    SolverType solverType;
    int timep;
};