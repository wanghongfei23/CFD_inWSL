
#pragma once
#include "block.hpp"
#include "oneDBnd.hpp"
#include "info.hpp"
#include <boost/circular_buffer.hpp>
#include "eigenSystem.hpp"




class SpaceDis
{
    public:
    SpaceDis(int n_,Data* data_,Data* rhs_
            ,std::shared_ptr<OneDBnd> bndL_,std::shared_ptr<OneDBnd> bndR_,Info* info);
    SpaceDis();
    void difference();
    void setOffset(int,int);
    void setMethod(EquationType ,DiffMethod);
    void setIDim(int);

    void setConstNorm(std::array<real,3>&&);

    //计时相关
    long long timep=0;
    


    private:
    
    Data* data;
    Data* rhs;
    Info* info;
    std::shared_ptr<OneDBnd> bndL,bndR;

    void calFlux();
    void calFluxConv(int);
    void calFluxBurgers(int);
    void calFluxEuler1D(int);
    void calFluxEuler2D(int);
    void calFluxAccuracyTest(int i);
    void (SpaceDis::*calTypeFlux)(int);

    real& at(int,int);
    real& fluxAt(int,int);
    std::vector<real> (SpaceDis::*reconLMethod)(int i);
    std::vector<real> (SpaceDis::*reconRMethod)(int i);
    std::vector<real> reconL(int);
    std::vector<real> reconLprim(int);
    std::vector<real> reconLChar1D(int i);
    std::vector<real> reconLChar2D(int i);
    std::vector<real> reconR(int);
    std::vector<real> reconRprim(int);
    std::vector<real> reconRChar1D(int i);
    std::vector<real> reconRChar2D(int i);


    std::vector<real> recon1DFaceCenter(int i);
    std::vector<real> recon2DFaceCenter(int i);


    real (*inter5) (std::array<real,5>)=nullptr;
    real (*inter5Positive) (std::array<real,5>)=nullptr;


    void difHCS();
    void difTraditional6();
    void dif2Order();
    void difMND6();
    void (SpaceDis::*difMethod)();
    
    int n,nVar,nPrim;
    int idim;
    int i0=0,offset=1;
    int center=0,centerOffset=1;
    int nHalf;

    std::array<real,3> norm;

    std::shared_ptr<Data> flux_d;
    std::shared_ptr<OneDBnd> fBndL,fBndR;
    EquationType fluxType;
    DiffMethod diffMethod=TRAD2;
    InterMethod interMethod=FIRSTORDER;


    //circular_buffer method
    //尝试用循环数组减小内存开销，觉得没啥意思算了
    //boost::circular_buffer<real> flux,var,fluxNode;
    // std::vector<boost::circular_buffer<real>> flux,var,fluxNode;
    // void initBuffer();
    // std::vector<real> getEulerPrim(int i);
    // std::vector<real> (SpaceDis::*variablesGetter)(int i)=&getEulerPrim;
    

};