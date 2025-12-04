// 包含头文件，用于访问Info类的声明和其他依赖项
#include "info.hpp"

// 类成员函数的实现

// 返回ghost cell的数量，这取决于使用的空间方法(spMethod)和差分方法(diffMethod)
int Info::nGhostCell()
{
    // 条件判断语句，根据不同的方法组合返回相应的ghost cell数量
    if(WCNS5==spMethod && TRAD6 == diffMethod) return 5;
    else if(WCNS5==spMethod && HDS6 == diffMethod) return 3;
    else if(WCNS5==spMethod && MND6 == diffMethod) return 4;
    else if(MUSCL==spMethod && TRAD2 == diffMethod) return 2;
    else if(FIRSTORDER==spMethod && TRAD2 == diffMethod) return 1;

    else
    {
        // 错误输出，当方法组合未定义时提示
        std::cout<<"Info error: undifined spMethod and diffMethod combination\n";
    }
    return 0;
}

// 返回通量点的数量，取决于使用的差分方法
int Info::nFluxPoint()
{
    // switch语句，根据不同差分方法返回相应的通量点数量
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

// 返回基本变量的数量，取决于控制方程类型和维度
int Info::nPrim()
{
    // switch语句，根据不同方程类型返回相应的基本变量数量
    switch (eqType)
    {
    case ACCURACYTEST:
    case LINEARCONV1D:
    case BURGERS1D:
        return 1;
        break;
    case EULER:
        {
            // 根据空间维度返回不同数量的变量
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

// 返回守恒变量的数量，取决于控制方程类型和维度
int Info::nCons()
{
    // switch语句，根据不同方程类型返回相应的守恒变量数量
    switch (eqType)
    {
    case ACCURACYTEST:
    case LINEARCONV1D:
    case BURGERS1D:
        return 1;
        break;
    case EULER:
        {
            // 根据空间维度返回不同数量的变量
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

// 计算并返回问题的空间维度
int Info::getDim()
{
    // 循环遍历三个空间方向，计算实际使用的维度数
    int dim=0;
    for (int i = 0; i < 3; i++)
    {
        // 如果某个方向的网格点数大于2，则该方向被使用
        if (iMax[i]>2) dim++;
    }
    return dim;
}

// 计算内部网格点的最大索引
std::array<int,3> Info::icMax()
{
    // 创建并初始化数组
    std::array<int,3> icMax;
    // 循环处理每个空间方向
    for (int i = 0; i < 3; i++)
    {
        // 根据网格点数设置内部点的最大索引
        if (iMax[i]<2) icMax[i]=1;
        else
        icMax[i]=iMax[i]-1;
    }
    return icMax;
}

// 返回默认的边界类型，取决于控制方程类型
BndType Info::defaultBndType()
{
    // switch语句，根据不同方程类型返回相应的默认边界类型
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

// 定义静态映射表，将方程类型映射到字符串
// static std::map<EquationType,std::string> fluxStr= {{LINEARCONV1D,"LINEARCONV1D"}
//                                         ,{BURGERS1D,"BURGERS1D"}
//                                         ,{EULER,"EULER"}
//                                         ,{ACCURACYTEST,"ACCURACYTEST"}};

// 【王鸿飞】begin输出命名系统cgns

// 定义静态映射表，将插值方法映射到字符串
static std::map<InterMethod,std::string> disStr={
    {FIRSTORDER,"FIRSTORDER"},
    {MUSCL,"MUSCL"},
    {WCNS5,"WCNS5"},
    {WCNSZ5,"WCNSZ5"},
    {WCNS5Char,"WCNS5Char"},
    {WCNSZ5Char,"WCNSZ5Char"},
    {WCNS5CONG,"WCNS5CONG"},
    {TCNSCongA,"TCNSCongA"},
    {WCNS5CONGZ,"WCNS5CONGZ"},
    {WCNS5CONGZCT4,"WCNS5CONGZCT4"},
    {WCNS5CONGZCT7,"WCNS5CONGZCT7"},
    {TCNS5,"TCNS5"},
    {TCNS5CT4,"TCNS5CT4"},
    {TCNS5CT7,"TCNS5CT7"},
    {LINEAR5,"LINEAR5"},
    {MUCSLIN5,"MUCSLIN5"},
    {INTERMAX,"INTERMAX"},
    //
    {whfTCNSN,"whfTCNSN"},
    {whfTCNSNA,"whfTCNSNA"},
    {whfTCNSNAS,"whfTCNSNAS"},
    {whfTCNSNS,"whfTCNSNS"},
    {whfTCNSNLAD,"whfTCNSNLAD"},
    //
    {congTCNS5CT5,"congTCNS5CT5"},
    {congTCNS5CT10,"congTCNS5CT10"},
    {whfAITCNSNS,"whfAITCNSNS"},
    {whfAITCNSNA,"whfAITCNSNA"},
    {whfAITCNSNAS_1,"whfAITCNSNAS_1"},
    {whfAITCNSNAS_2,"whfAITCNSNAS_2"},
    {whfAITCNSNLADS,"whfAITCNSNLADS"},
    {whfAITCNSNAZS,"whfAITCNSNAZS"},
    {whfAITCNSNmyASF002_1,"whfAITCNSNmyASF002_1"},
    {whfAITCNSNmyASF002_2,"whfAITCNSNmyASF002_2"},
    {whfAITCNSNmyASF002_ai1,"whfAITCNSNmyASF002_ai1"},
    {whfzycTCNSNmyASF002_1,"whfzycTCNSNmyASF002_1"},
    {whfCOMPARE,"whfCOMPARE"},
    {new_TCNS5,"TENO"},
    {new_WHFTCNSA,"TENO-A"},
    {new_WCNS5CONGZ,"TENO-S"},
    {new_WHFTCNSASF203_NoS,"TENO-AS-myF203_NoS"},
    {new_WHFTCNSASF202_NoS,"TENO-AS-myF202_NoS"}
};
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
    {0,"2D_Riemann"},
    {1,"2D_Riemann_vortex"},
    {2,"implosion"},
    {3,"RTI"},
    {4,"Double_Mach"},
    {5,"2D_Riemann_another"},
    {6,"KHI"}
};


// 生成输出文件名   // 原命名规则
// std::string Info::filename()
// {
//     // 字符串拼接，包含方程类型、插值方法和当前时间
//     return fluxStr[eqType]+disStr[spMethod]+std::format("t={:.4f}.cgns",t);
// }

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

// 【王鸿飞】end输出命名系统cgns

// 获取守恒变量名称列表
std::vector<std::string> Info::getVarNameListCons()
{
    // 创建字符串向量
    std::vector<std::string> res;
    // 预留空间以提高性能
    res.reserve(nCons());
    // 根据维度和方程类型添加不同的变量名
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

// 获取基本变量名称列表
std::vector<std::string> Info::getVarNameListPrim()
{
    // 创建字符串向量
    std::vector<std::string> res;
    // 预留空间以提高性能
    res.reserve(nPrim());
    // 根据维度和方程类型添加不同的变量名
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

// 获取RHS变量名称列表
std::vector<std::string> Info::getVarNameListRhs()
{
    // 创建字符串向量
    std::vector<std::string> res;
    // 预留空间以提高性能
    res.reserve(nCons());
    // 根据维度和方程类型添加不同的变量名
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

// 计算指定方向上的网格步长
real Info::geth(int idim)
{
    // 根据是否使用恒定步长返回相应的值
    if(constH) return interval;
    // 计算实际的网格步长
    return (calZone[2*idim+1]-calZone[2*idim])/(iMax[idim]-1);
}

// Info类的构造函数，初始化对象的维度信息
Info::Info()
{
    // 调用getDim()函数计算并设置空间维度
    dim=getDim();
}