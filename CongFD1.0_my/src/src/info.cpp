/**
 * @file info.cpp
 * @brief Info类的实现文件
 */

#include "info.hpp"

/**
 * @brief 获取虚拟点单元数
 * @return 虚拟点单元数
 */
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

/**
 * @brief 获取通量点数
 * @return 通量点数
 */
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

/**
 * @brief 获取原始变量数
 * @return 原始变量数
 */
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

/**
 * @brief 获取守恒变量数
 * @return 守恒变量数
 */
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


/**
 * @brief 获取问题维度
 * @return 问题维度
 */
int Info::getDim()
{
    int dim=0;
    for (int i = 0; i < 3; i++)
    {
        if (iMax[i]>2) dim++;
    }
    return dim;
}

/**
 * @brief 获取内部网格最大索引
 * @return 内部网格最大索引数组
 */
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


/**
 * @brief 获取默认边界类型
 * @return 默认边界类型
 */
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
    // 【王鸿飞】begin-1命名
    {WCNS5,"WENO-JS"},
    {WCNSZ5,"WENO-Z"},
    {TCNS5,"TENO"},
    {WCNS5CONGZ,"TENO-S"},
    {WHFTCNSA,"TENO-myA"},
    {WHFTCNSASF002,"TENO-AS-myF002"},
    {WHFTCNSAH002,"TENO-AS-myH002"},
    {WHFTCNSASF102,"TENO-AS-myF102"},
    {WHFTCNSASF103,"TENO-AS-myF103"},
    {WHFTCNSASF102_reciprocal,"TENO-AS-myF102_reciprocal"},
    {WHFTCNSASF103_reciprocal,"TENO-AS-myF103_reciprocal"},
    {WHFTCNSAS_fx,"TENO-AS-fx"},
    {WHFTCNSAS_initial,"TENO-AS-initial"},
    {WHFTCNSAS_approx_1,"TENO-AS-approx_1"},
    {WHFTCNSAS_fx_real,"TENO-AS-fx-real"},
    {WHFTCNSAS_approx_2,"TENO-AS-approx_2"},
    {WHFTCNSASF202_2S,"TENO-AS-myF202_2S"},
    {WHFTCNSASF202_NoS,"TENO-AS-myF202_NoS"},
    {WHFTCNSASF203_NoS,"TENO-AS-myF203_NoS"},
    {temp011,"temp_name_011"},
    {temp012,"temp_name_012"},
    {temp013,"temp_name_013"},
    {temp014,"temp_name_014"},
    {temp015,"temp_name_015"},
    {temp016,"temp_name_016"},
    {temp017,"temp_name_017"},
    {temp018,"temp_name_018"},
    {temp019,"temp_name_019"}
    // 【王鸿飞】end-1命名
};

// 新命名系统

// 定义静态映射表，将算例映射到字符串
static std::map<int,std::string> exampleStr1D={
    {0,"Sod"},
    {1,"ShuOsher"},
    {2,"Lax"},
    {3,"sedov"},
    {4,"Woodward_Colella"},
    {5,"Double_sparse_wave"}
};

static std::map<int,std::string> exampleStr2D={
    {0,"2D_Riemann_1"},
    {1,"2D_Riemann_2"},
    {2,"implosion"},
    {3,"RTI"},
    {4,"Double_Mach"},
    {5,"2D_Riemann_3"},
    {6,"KHI"}
};

/**
 * @brief 生成文件名
 * @return 生成的文件名
 */
std::string Info::filename()
{
    // 字符串拼接，包含空间离散方法、插值方法、算例信息和当前时间
    std::string caseName = "unknown";
    
    if (dim == 1) {
        auto it = exampleStr1D.find(nCase);
        if (it != exampleStr1D.end()) {
            caseName = it->second;
        }
    } else if (dim == 2) {
        auto it = exampleStr2D.find(nCase);
        if (it != exampleStr2D.end()) {
            caseName = it->second;
        }
    }
    
    // 添加网格信息到文件名
    std::string gridInfo = "";
    if (dim == 1) {
        gridInfo = std::to_string(iMax[0]);
    } else if (dim == 2) {
        gridInfo = std::to_string(iMax[0]) + "x" + std::to_string(iMax[1]);
    } else if (dim == 3) {
        gridInfo = std::to_string(iMax[0]) + "x" + std::to_string(iMax[1]) + "x" + std::to_string(iMax[2]);
    }
    
    return caseName + " - " + disStr[interMethod] + " - " + gridInfo + " - " + std::format("t={:.4f}.cgns",t);
}


// 原命名系统
// std::string Info::filename()
// {
//     return fluxStr[eqType]+disStr[spMethod]+std::format("t={:.4f}.cgns",t);
// }


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