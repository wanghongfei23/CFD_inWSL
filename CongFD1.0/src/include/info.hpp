#pragma once
#include "macro.hpp"

class Info
{
    public:
    InterMethod spMethod=WCNS5;
    EquationType eqType=EULER;
    int nCase=2;
    real t=0;
    real dt=0.0001;
    int step=0;
    int endStep=25000;
    int dim=1;


    int outputInterval=100;
    real outputDt=0.01;
    real outputT=0;

    real CFL=0.5;
    bool fixedtimeSteps=true;

    //for implicit solver
    real implicitCFL=0.01;
    real maxImplicitStep=100;



    DiffMethod diffMethod=TRAD6;
    InterMethod interMethod=WCNS5;//只影响插值权
    SourceType sourceType=SOURCENULL;

    std::array<int,3> iMax{201,201,2};
    std::array<double,6> calZone{0,0.3,0,0.3,0,2};

    Info();
    int nGhostCell();
    int nFluxPoint();
    int nPrim();
    int nCons();
    int getDim();

    bool constH=true;
    real interval=0;
    real geth(int);


    BndType defaultBndType();
    std::array<int,3> icMax();
    std::string filename();

    std::vector<std::string> getVarNameListCons();
    std::vector<std::string> getVarNameListPrim();
    std::vector<std::string> getVarNameListRhs();
};