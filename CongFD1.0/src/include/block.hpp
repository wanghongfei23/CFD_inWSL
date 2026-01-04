#pragma once
#include"data.hpp"
#include <cgnslib.h>

class Block
{
    public:
    void outputCgns();
    real operator()(int,int);
    std::array<int,3> getICMax();
    std::array<int,3> getIMax();
    real getMinDh(int i);
    int getDim();
    std::vector<real> getCellCoor(int idim);
    std::vector<real> getVertexCoor(int idim);
    std::vector<real> getCellInterval(int);
    


    protected:
    friend class Initializer;
    Data coorVer;
    Data coorCel;
    Data intervalCel;
    int nVer,nCel;
    bool inited=false;
    int dim;
    std::array<int,3> iMax,icMax;

};