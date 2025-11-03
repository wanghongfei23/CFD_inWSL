#pragma once

#include"data.hpp"

class Equation
{
    public:
    Equation(){};
    void consToPrim();
    
    Data* getCons();
    Data* getPrim();
    Data* getRhs();

    protected:
    friend class Initializer;

    void consToPrimEuler1D();
    void consToPrimEuler2D();
    bool inited=false;
    int n,nCons,nPrim;
    int dim;
    Data* cons,*rhs,*prim;
    EquationType type;
};