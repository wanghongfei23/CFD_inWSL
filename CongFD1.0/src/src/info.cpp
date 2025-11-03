#include "info.hpp"

int Info::nGhostCell()
{
    
    if(WCNS5==spMethod && TRAD6 == diffMethod) return 5;
    else if(WCNS5==spMethod && HDS6 == diffMethod) return 3;
    else if(WCNS5==spMethod && MND6 == diffMethod) return 4;
    else if(MUSCL==spMethod && TRAD2 == diffMethod) return 2;
    else if(FIRSTORDER==spMethod && TRAD2 == diffMethod) return 1;

    else
    {
        std::cout<<"Info error: undifined spMethod and diffMethod combination\n";
    }
    return 0;
}
int Info::nFluxPoint()
{
    switch (diffMethod)
    {
    case TRAD2:
    case HDS6:
        return 0;
        break;
    case MND6:
        return 1;
        break;
    case TRAD6:
        return 2;
        break;
    
    default:
        return 0;
        break;
    }
}

int Info::nPrim()
{
    switch (eqType)
    {
    case ACCURACYTEST:
    case LINEARCONV1D:
    case BURGERS1D:
        return 1;
        break;
    case EULER:
        {
            if(dim==1) return 3;
            if(dim==2) return 4;
            else return 0;
        }
        break;
    
    default:
        std::cout<<"Info error: undifined eqType in nPrim()\n";
        return 0;
        
        break;
    }
}
int Info::nCons()
{
    switch (eqType)
    {
    case ACCURACYTEST:
    case LINEARCONV1D:
    case BURGERS1D:
        return 1;
        break;
    case EULER:
        {
            if(dim==1) return 3;
            if(dim==2) return 4;
            else return 0;
        }
        break;
    
    default:
        std::cout<<"Info error: undifined eqType in nCons()\n";
        return 0;
        
        break;
    }
}


int Info::getDim()
{
    int dim=0;
    for (int i = 0; i < 3; i++)
    {
        if (iMax[i]>2) dim++;
    }
    return dim;
}

std::array<int,3> Info::icMax()
{
    std::array<int,3> icMax;
    for (int i = 0; i < 3; i++)
    {
        if (iMax[i]<2) icMax[i]=1;
        else
        icMax[i]=iMax[i]-1;
    }
    return icMax;
}


BndType Info::defaultBndType()
{
    switch (eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        return PERIODIC1D;
        break;
    
    default:
        return TYPENULL;
        break;
    }
}

static std::map<EquationType,std::string> fluxStr= {{LINEARCONV1D,"LINEARCONV1D"}
                                        ,{BURGERS1D,"BURGERS1D"}
                                        ,{EULER,"EULER"}
                                        ,{ACCURACYTEST,"ACCURACYTEST"}};
static std::map<InterMethod,std::string> disStr={
    {FIRSTORDER,"FIRSTORDER"},
    {MUSCL,"MUSCL"},
    {WCNS5,"WCNS5"},
    {WCNSZ5,"WCNSZ5"}
};

std::string Info::filename()
{
    return fluxStr[eqType]+disStr[spMethod]+std::format("t={:.4f}.cgns",t);
}

std::vector<std::string> Info::getVarNameListCons()
{
    std::vector<std::string> res;
    res.reserve(nCons());
    if(dim==1)
    {
        if(eqType==EULER)
        {
            res.push_back("rho");
            res.push_back("rhoU");
            res.push_back("rhoE");
        }
        else res.push_back("u");
    }
    else if (dim==2)
    {
        if(eqType==EULER)
        {
            res.push_back("rho");
            res.push_back("rhoU");
            res.push_back("rhoV");
            res.push_back("rhoE");
        }
    }
    return res;
}
std::vector<std::string> Info::getVarNameListPrim()
{
    std::vector<std::string> res;
    res.reserve(nPrim());
    if(dim==1)
    {
        if(eqType==EULER)
        {
            res.push_back("Density");
            res.push_back("XVelocity");
            res.push_back("Pressure");
            // res.push_back("Density");
            // res.push_back("u+c");
            // res.push_back("u-c");
        }
        else res.push_back("u");
    }
    else if (dim==2)
    {
        if(eqType==EULER)
        {
            res.push_back("Density");
            res.push_back("XVelocity");
            res.push_back("YVelocity");
            res.push_back("Pressure");
        }
    }
    return res;
}

std::vector<std::string> Info::getVarNameListRhs()
{
    std::vector<std::string> res;
    res.reserve(nCons());
    if(dim==1)
    {
        if(eqType==EULER)
        {
            res.push_back("RHS-rho");
            res.push_back("RHS-rhoU");
            res.push_back("RHS-rhoE");
        }
        else res.push_back("RHS-u");
    }
    else if (dim==2)
    {
        if(eqType==EULER)
        {
            res.push_back("RHS-rho");
            res.push_back("RHS-rhoU");
            res.push_back("RHS-rhoV");
            res.push_back("RHS-rhoE");
        }
    }
    return res;
}


real Info::geth(int idim)
{
    if(constH) return interval;
    return (calZone[2*idim+1]-calZone[2*idim])/(iMax[idim]-1);
}

Info::Info()
{
    dim=getDim();
}