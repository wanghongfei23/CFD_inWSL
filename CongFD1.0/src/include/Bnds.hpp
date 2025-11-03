#pragma once
#include "info.hpp"
#include "oneDBnd.hpp"
#include "data.hpp"





class Bnds
{
    public:
    std::array<std::shared_ptr<OneDBnd>,2> getOneDBnd(int,int,int);
    void update();

    private:
    friend class Initializer;
    int dim;
    std::array<int,3> iMax;
    //idim*
    std::vector<std::shared_ptr<OneDBnd>> oneDBnds;//JK*2+IK*2+IJ*2
};