#pragma once

#include "data.hpp"
#include "info.hpp"


class fluxSolver
{
    public:
    fluxSolver(std::shared_ptr<Data> valsR_,std::shared_ptr<Data> fluxes,Info* info_,int idim_);
    virtual void fluxSolve()=0;

    protected:
    std::shared_ptr<Data> valsR;
    std::shared_ptr<Data> fluxes;
    Info* info;
    int idim;


};