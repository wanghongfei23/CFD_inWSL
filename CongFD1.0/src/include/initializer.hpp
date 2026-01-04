#pragma once
#include"block.hpp"
#include"info.hpp"
#include"equation.hpp"
#include"Bnds.hpp"
#include"sp_distributor.hpp"

class Initializer
{
    public:
    Initializer();
    Initializer(Info*);
    void solInit(Block*,Data*);
    void initUniformBlock(Block*);
    void initEqution(Equation*,Block*);



    void initBnds(Bnds* bnds,Equation*,std::array<int,3> iMax,Block* block);
    void initDoubleMachBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block);


    void initSpDistributor(SpDistributor*,Equation*,Block*,Bnds*);


    private:
    Info* info;
    
    
};