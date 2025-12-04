#include "solverType.hpp"

/**
 * @brief 拷贝构造函数实现
 * @param origin 被拷贝的SolverType对象
 */
SolverType::SolverType(const SolverType& origin)
{
    info = origin.info;
    reconer = origin.reconer;
    solvePointSolver = origin.solvePointSolver;
    fluxPointSolver = origin.fluxPointSolver;
    differ = origin.differ;
}

/**
 * @brief 带参数的构造函数实现
 * @param info_ 指向Info对象的指针，包含求解器配置信息
 */
SolverType::SolverType(Info* info_)
{
    info = info_;
    initReconer();
    initFluxPointSolver();
    initSolvePointSolver();
    initDiffer();
}

/**
 * @brief 初始化重构器的实现
 * 
 * 根据方程类型和插值方法选择合适的重构器进行初始化
 */
void SolverType::initReconer()
{
    switch (info->eqType) {
    case LINEARCONV1D:
    case BURGERS1D:
        switch (info->interMethod) {
        case WCNS5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<weno5_JSchen>());
            break;
        case WCNSZ5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<weno5_Z>());
            break;
        case TCNS5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_Z>());
            break;
        case TCNS5CT4:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_ZCT4>());
            break;
        case TCNS5CT7:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_ZCT7>());
            break;
        case WCNS5CONGZ:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_CongZ>());
            break;
        case WCNS5CONGZCT4:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_CongZCT4>());
            break;
        case WCNS5CONGZCT7:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<Teno5_CongAA>());
            break;

        //【王鸿飞】begin插值格式开发（1.0）手搓
        case whfTCNSN:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_TcnsN>());
            break;

        case whfTCNSNA:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_TcnsN_A>());
            break;

        case whfTCNSNAS:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_TcnsN_AS>());
            break;

        case whfTCNSNS:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_TcnsN_S>());
            break;

        case whfTCNSNLAD:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_TcnsN_LAD>());
            break;
        //【王鸿飞】end插值格式开发（1.0）手搓

        //【王鸿飞】begin插值格式开发（2.0）AI-原格式

        case congTCNS5CT5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<congTcns5_ZCT5>());
            break;

        case congTCNS5CT10:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<congTcns5_ZCT10>());
            break;
        
        case whfAITCNSNS:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_S_1>());
            break;
            
        case whfAITCNSNA:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_A_2>());
            break;

        case whfAITCNSNAS_1:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_AS_1>());
            break;

        case whfAITCNSNAS_2:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_AS_2>());
            break;

        case whfAITCNSNLADS:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_LADS_yuanbao>());
            break;

        case whfAITCNSNAZS:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_AZS_yuanbao>());
            break;

        //【王鸿飞】end插值格式开发（2.0）AI-原格式

        // 【王鸿飞】begin插值格式开发（3.0）AI-新格式

        case whfAITCNSNmyASF002_1:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_myASF002_1>());
            break;

        case whfAITCNSNmyASF002_2:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_myASF002_2>());
            break;
            
        case whfAITCNSNmyASF002_ai1:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_ai_TcnsN_myASF002_ai1>());
            break;



        // 【王鸿飞】end插值格式开发（3.0）AI-新格式

        // 【王鸿飞】begin插值格式（4.0）手搓及优化

        case whfzycTCNSNmyASF002_1:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_zyc_TcnsN_myASF002_1>());
            break;

        case whfCOMPARE:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_TcnsN_COMPARE>());
            break;

        case new_TCNS5:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<new_Teno5_Z>());
            break;


        case new_WHFTCNSA:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<new_whf_TCNS_A>());
            break;


        case new_WCNS5CONGZ:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<new_Teno5_CongZ>());
            break;


        case new_WHFTCNSASF203_NoS:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_TCNS_AS_myF203_NoS>());
            break;

        case new_WHFTCNSASF202_NoS:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<whf_TCNS_AS_myF202_NoS>());
            break;




        // 【王鸿飞】end插值格式（4.0）手搓及优化


        default:
            reconer = pro::make_proxy<ProxyDataManipulator>(
                Recon5OrderFaceCenter<weno5_JSchen>());
            break;
        }
        break;
    case EULER:
        if (info->dim == 1) {
            switch (info->interMethod) {
            case WCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<weno5_JSchen>());
                break;
            case WCNSZ5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<weno5_Z>());
                break;
            case TCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_Z>());
                break;
            case TCNS5CT4:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_ZCT4>());
                break;
            case TCNS5CT7:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_ZCT7>());
                break;
            case WCNS5CONGZ:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_CongZ>());
                break;
            case WCNS5CONGZCT4:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_CongZCT4>());
                break;
            case WCNS5CONGZCT7:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<Teno5_CongAA>());
                break;
            case MUCSLIN5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<musclIn5>());
                    break;

            //【王鸿飞】begin插值格式开发（1.0）手搓
            case whfTCNSN:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_TcnsN>());
                break;

            case whfTCNSNA:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_TcnsN_A>());
                break;

            case whfTCNSNAS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_TcnsN_AS>());
                break;

            case whfTCNSNS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_TcnsN_S>());
                break;

            case whfTCNSNLAD:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_TcnsN_LAD>());
                break;
            //【王鸿飞】end插值格式开发（1.0）手搓

            //【王鸿飞】begin插值格式开发（2.0）AI-原格式

            case congTCNS5CT5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<congTcns5_ZCT5>());
                break;

            case congTCNS5CT10:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<congTcns5_ZCT10>());
                break;
                
            case whfAITCNSNS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_S_1>());
                break;
                
            case whfAITCNSNA:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_A_2>());
                break;
                
            case whfAITCNSNAS_1:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_AS_1>());
                break;
                
            case whfAITCNSNAS_2:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_AS_2>());
                break;
                
            case whfAITCNSNLADS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_LADS_yuanbao>());
                break;

            case whfAITCNSNAZS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_AZS_yuanbao>());
                break;



            //【王鸿飞】end插值格式开发（2.0）AI-原格式

            // 【王鸿飞】begin插值格式开发（3.0）AI-新格式

            case whfAITCNSNmyASF002_1:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_myASF002_1>());
                break;

            case whfAITCNSNmyASF002_2:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_myASF002_2>());
                break;
                
            case whfAITCNSNmyASF002_ai1:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_ai_TcnsN_myASF002_ai1>());
                break;


            // 【王鸿飞】end插值格式开发（3.0）AI-新格式
            

            // 【王鸿飞】begin插值格式（4.0）手搓及优化

            case whfzycTCNSNmyASF002_1:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_zyc_TcnsN_myASF002_1>());
                break;

            case whfCOMPARE:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_TcnsN_COMPARE>());
                break;


            case new_TCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<new_Teno5_Z>());
                break;


            case new_WHFTCNSA:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<new_whf_TCNS_A>());
                break;


            case new_WCNS5CONGZ:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<new_Teno5_CongZ>());
                break;


            case new_WHFTCNSASF203_NoS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_TCNS_AS_myF203_NoS>());
                break;

            case new_WHFTCNSASF202_NoS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<whf_TCNS_AS_myF202_NoS>());
                break;





            // 【王鸿飞】end插值格式（4.0）手搓及优化



            default:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order1DEulerEig<weno5_JSchen>());
                break;
            }
        } else if (info->dim == 2) {
            switch (info->interMethod) {
            case WCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<weno5_JSchen>());
                break;
            case WCNSZ5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<weno5_Z>());
                break;
            case TCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_Z>());
                break;
            case TCNS5CT4:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_ZCT4>());
                break;
            case TCNS5CT7:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_ZCT7>());
                break;
            case WCNS5CONGZ:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_CongZ>());
                break;
            case WCNS5CONGZCT4:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_CongZCT4>());
                break;
            case WCNS5CONGZCT7:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<Teno5_CongZCT7>());
                break;
            case MUCSLIN5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<musclIn5>());
                break;
                
            //【王鸿飞】begin插值格式开发（1.0）手搓
            case whfTCNSN:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_TcnsN>());
                break;

            case whfTCNSNA:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_TcnsN_A>());
                break;

            case whfTCNSNAS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_TcnsN_AS>());
                break;

            case whfTCNSNS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_TcnsN_S>());
                break;

            case whfTCNSNLAD:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_TcnsN_LAD>());
                break;
            //【王鸿飞】end插值格式开发（1.0）手搓

            //【王鸿飞】begin插值格式开发（2.0）AI-原格式

            case congTCNS5CT5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<congTcns5_ZCT5>());
                break;

            case congTCNS5CT10:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<congTcns5_ZCT10>());
                break;
                
            case whfAITCNSNS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_S_1>());
                break;
                
            case whfAITCNSNA:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_A_2>());
                break;
                
            case whfAITCNSNAS_1:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_AS_1>());
                break;

            case whfAITCNSNAS_2:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_AS_2>());
                break;
                
            case whfAITCNSNLADS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_LADS_yuanbao>());
                break;
                
            case whfAITCNSNAZS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_AZS_yuanbao>());
                break;


            //【王鸿飞】end插值格式开发（2.0）AI-原格式


            // 【王鸿飞】begin插值格式开发（3.0）AI-新格式

            case whfAITCNSNmyASF002_1:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_myASF002_1>());
                break;

            case whfAITCNSNmyASF002_2:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_myASF002_2>());
                break;
                
            case whfAITCNSNmyASF002_ai1:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_ai_TcnsN_myASF002_ai1>());
                break;

            // 【王鸿飞】end插值格式开发（3.0）AI-新格式




            // 【王鸿飞】begin插值格式（4.0）手搓及优化

            case whfzycTCNSNmyASF002_1:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_zyc_TcnsN_myASF002_1>());
                break;

            case whfCOMPARE:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_TcnsN_COMPARE>());
                break;


            case new_TCNS5:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<new_Teno5_Z>());
                break;


            case new_WHFTCNSA:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<new_whf_TCNS_A>());
                break;


            case new_WCNS5CONGZ:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<new_Teno5_CongZ>());
                break;


            case new_WHFTCNSASF203_NoS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_TCNS_AS_myF203_NoS>());
                break;

            case new_WHFTCNSASF202_NoS:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<whf_TCNS_AS_myF202_NoS>());
                break;





            // 【王鸿飞】end插值格式（4.0）手搓及优化



            default:
                reconer = pro::make_proxy<ProxyDataManipulator>(
                    Recon5Order2DEulerEig<weno5_JSchen>());
                break;
            }
        }
        break;
    default:
        break;
    }
}

void SolverType::initSolvePointSolver()
{
    switch (info->eqType) {
    case LINEARCONV1D:
        solvePointSolver = pro::make_proxy<ProxySolvePointSolver>(
            SolvePointSolver<1, fluxSolveLinearConv>());
        break;
    case BURGERS1D:
        solvePointSolver = pro::make_proxy<ProxySolvePointSolver>(
            SolvePointSolver<1, fluxSolveBurgers>());
        break;
    case EULER:
        if (info->dim == 1) {
            solvePointSolver = pro::make_proxy<ProxySolvePointSolver>(
                SolvePointSolver<3, fluxSolveEuler1D>());
        } else if (info->dim == 2) {
            solvePointSolver = pro::make_proxy<ProxySolvePointSolver>(
                SolvePointSolver<4, fluxSolveEuler2D>());
        }
        break;
    default:
        break;
    }
}

void SolverType::initFluxPointSolver()
{
    switch (info->eqType) {
    case EULER:
        switch (info->fluxMethod) {
        case ROE:
            if (info->dim == 1) {
                fluxPointSolver = pro::make_proxy<ProxyFluxPointSolver>(
                    FluxPointSolver<3, roeFlux1D2>());
            } else if (info->dim == 2) {
                fluxPointSolver = pro::make_proxy<ProxyFluxPointSolver>(
                    FluxPointSolver<4, roeFlux2DSym>());
            }
            break;
        case HLLC:
            if (info->dim == 1) {
                fluxPointSolver = pro::make_proxy<ProxyFluxPointSolver>(
                    FluxPointSolver<3, HLLCFlux1D>());
            } else if (info->dim == 2) {
                fluxPointSolver = pro::make_proxy<ProxyFluxPointSolver>(
                    FluxPointSolver<4, HLLCFlux2D2>());
            }
            break;

        default:
            break;
        }
        break;
    default:
        break;
    }
}

void SolverType::initDiffer()
{
    switch (info->eqType) {
    case EULER:
        switch (info->diffMethod) {
        case TRAD2:
            differ = pro::make_proxy<ProxyDiffer>(Differ());
            solvePointSolver->setNlNr(0, 0);
            break;
        case TRAD6:
            differ = pro::make_proxy<ProxyDiffer>(ExplicitDif6());
            solvePointSolver->setNlNr(0, 0);
            break;
        case MND6:
            differ = pro::make_proxy<ProxyDiffer>(MidNodeAndNodeDif6());
            solvePointSolver->setNlNr(1, 1);
            break;
        default:
            break;
        }
        break;
    default:
        break;
    }
}