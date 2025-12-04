#pragma once
#include "macro.hpp"
#include "data.hpp"
#include "info.hpp"
#include "oneDBnd.hpp"

class Reconstructor
{
    public:
    Reconstructor(int n_,Data* vals_,std::shared_ptr<Data> valsR_,std::shared_ptr<OneDBnd> bndL_,std::shared_ptr<OneDBnd> bndR_,Info* info_,int i0_, int offset_);
    virtual void recon()=0;
    
    std::shared_ptr<Data> valsR;

    protected:
    Data* vals;
    real& at(int i,int ivar);
    int i0,offset;
    std::shared_ptr<OneDBnd> bndL,bndR;
    Info* info;
    int n,nval;
};

class ReconEigen1DEuler : public Reconstructor
{
    public:
    void recon() override;
};