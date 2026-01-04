#include "equation.hpp"

void Equation::consToPrim()
{
    switch (type)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        std::copy(cons->begin(),cons->end(),prim->begin());
        break;
    case EULER:
        {
            if (dim==1) consToPrimEuler1D();
            else if(dim==2) 
            {
                consToPrimEuler2D();
            }
        }
        break;
    
    default:
        break;
    }
}

void Equation::consToPrimEuler1D()
{
    if(nCons!=3,nPrim!=3)
    {
        std::cout<<"Equation error: Euler 1d equation variable number error \n";
    }

    for (int i = 0; i < n; i++)
    {
        real r=(*cons)(i,0);
        real ru=(*cons)(i,1);
        real rE=(*cons)(i,2);
        real u=ru/r;
        real E=rE/r;
        real e=-u*u/2+E;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=p;
    }
}

void Equation::consToPrimEuler2D()
{
    if(nCons!=3,nPrim!=4)
    {
        std::cout<<"Equation error: Euler 1d equation variable number error \n";
    }

    for (int i = 0; i < n; i++)
    {
        real r=(*cons)(i,0);
        real ru=(*cons)(i,1);
        real rv=(*cons)(i,2);
        real rE=(*cons)(i,3);
        real u=ru/r;
        real v=rv/r;
        real E=rE/r;
        real q2=(u*u+v*v)/2;
        real e=-q2+E;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=v;
        (*prim)(i,3)=p;
    }
}


Data* Equation::getPrim()
{
    return prim;
}

Data* Equation::getCons()
{
    return cons;
}

Data* Equation::getRhs()
{
    return rhs;
}