#include "SourceTerm.hpp"

void SourceTerm::calGravitySource()
{
    //here it is still only for R-T instability case
    if(nprim!=4) 
    {
        std::cout<<"SourceTerm error: nprim incorrect\n";
    }
    real r,u,v,p;
    for (int i = 0; i < n; i++)
    {
        r=(*prim)(i,0);
        u=(*prim)(i,1);
        v=(*prim)(i,2);
        p=(*prim)(i,3);
        (*rhs)(i,2)-=r;
        (*rhs)(i,3)-=r*v;
    }
    
}
void SourceTerm::nothingHappened()
{

}


SourceTerm::SourceTerm(Data* prim_,Data* rhs_,Info* info_)
{
    prim=prim_;
    rhs=rhs_;
    info=info_;

    auto iMax=info->icMax();
    n=iMax[0]*iMax[1]*iMax[2];
    nprim=info->nPrim();
    nCons=info->nCons();
    switch (info->sourceType)
    {
    case GRAVITY:
        calSourceMethod=(&SourceTerm::calGravitySource);
        break;
    
    default:
        calSourceMethod=(&SourceTerm::nothingHappened);
        break;
    }
}

void SourceTerm::calSource()
{
    (this->*calSourceMethod)();
}