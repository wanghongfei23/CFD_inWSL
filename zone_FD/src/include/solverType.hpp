#include <variant>

#include "differ.hpp"
#include "fluxPointFlux.hpp"
#include "fluxSchemes.hpp"
#include "proxy.h"
#include "reconstructor.hpp"
#include "reconstructor5order.hpp"
#include "solvePointFlux.hpp"

/// @brief 定义求解器的Solve方法调度
PRO_DEF_MEM_DISPATCH(Solve, solve);

/// @brief 定义求解器的Init方法调度
PRO_DEF_MEM_DISPATCH(Init, init);

/// @brief 定义设置法向量常量的方法调度
PRO_DEF_MEM_DISPATCH(SetConstNorm, setConstNorm);

/// @brief 定义获取数据的方法调度
PRO_DEF_MEM_DISPATCH(GetData, getData);

/// @brief 定义检查方法的调度
PRO_DEF_MEM_DISPATCH(Check, check);

/// @brief 定义设置NlNr的方法调度
PRO_DEF_MEM_DISPATCH(SetNlNr, setNlNr);

/// @brief 定义设置常量H的方法调度
PRO_DEF_MEM_DISPATCH(SetConstantH, setConstantH);

/**
 * @brief 数据操作代理结构体
 * 
 * 定义了数据操作代理的接口，包括求解、初始化、设置法向量等操作
 */
struct ProxyDataManipulator
    : pro::facade_builder ::add_convention<Solve, void()>::
          add_convention<Check, void()>::
              add_convention<Init, void(DataReader, std::shared_ptr<OneDBnd>, std::shared_ptr<OneDBnd>)>::
                  add_convention<SetConstNorm, void(std::array<real, 3>)>::
                      add_convention<GetData, std::shared_ptr<Data>()>::
                          support_copy<pro::constraint_level::nontrivial>::build { };

/**
 * @brief 求解点求解器代理结构体
 * 
 * 定义了求解点求解器代理的接口，扩展了ProxyDataManipulator的功能
 */
struct ProxySolvePointSolver
    : pro::facade_builder ::add_convention<Solve, void()>::
          add_convention<Check, void()>::
              add_convention<Init, void(DataReader, std::shared_ptr<OneDBnd>, std::shared_ptr<OneDBnd>)>::
                  add_convention<SetConstNorm, void(std::array<real, 3>)>::
                      add_convention<SetNlNr, void(int, int)>::
                          add_convention<GetData, std::shared_ptr<Data>()>::
                              support_copy<pro::constraint_level::nontrivial>::build { };

/**
 * @brief 差分代理结构体
 * 
 * 定义了差分代理的接口，用于处理数据差分操作
 */
struct ProxyDiffer
    : pro::facade_builder ::
          add_convention<Check, void()>::
              add_convention<Init, void(std::shared_ptr<Data>, std::shared_ptr<Data>, DataReader)>::
                  add_convention<SetConstantH, void(real)>::
                      add_convention<Solve, void()>::
                          support_copy<pro::constraint_level::nontrivial>::build { };

/**
 * @brief 通量点求解器代理结构体
 * 
 * 定义了通量点求解器代理的接口
 */
struct ProxyFluxPointSolver
    : pro::facade_builder ::
          add_convention<Check, void()>::
              add_convention<Solve, void()>::
                  add_convention<Init, void(std::shared_ptr<Data>, int, int)>::
                      add_convention<SetConstNorm, void(std::array<real, 3>)>::
                          add_convention<GetData, std::shared_ptr<Data>()>::
                              support_copy<pro::constraint_level::nontrivial>::build { };

/**
 * @brief 求解器类型类
 * 
 * SolverType类负责管理不同类型的求解器组件，包括重构器、求解点求解器、
 * 通量点求解器和差分器，并根据配置信息初始化这些组件。
 */
class SolverType {
public:
    /**
     * @brief 带参数的构造函数
     * @param info_ 指向Info对象的指针，包含求解器配置信息
     */
    SolverType(Info* info_);
    
    /**
     * @brief 默认构造函数
     */
    SolverType() {};
    
    /**
     * @brief 拷贝构造函数
     * @param origin 被拷贝的SolverType对象
     */
    SolverType(const SolverType& origin);
    
    pro::proxy<ProxyDataManipulator> reconer;          ///< 重构器代理
    pro::proxy<ProxySolvePointSolver> solvePointSolver; ///< 求解点求解器代理
    pro::proxy<ProxyFluxPointSolver> fluxPointSolver;   ///< 通量点求解器代理
    pro::proxy<ProxyDiffer> differ;                     ///< 差分器代理

private:
    Info* info; ///< 指向Info对象的指针，包含求解器配置信息
    
    /**
     * @brief 初始化重构器
     */
    void initReconer();
    
    /**
     * @brief 初始化求解点求解器
     */
    void initSolvePointSolver();
    
    /**
     * @brief 初始化通量点求解器
     */
    void initFluxPointSolver();
    
    /**
     * @brief 初始化差分器
     */
    void initDiffer();
};

/*世界上最丑陋的代码*/
// using ReconstructorType = std::variant<
//     Recon1Order, Recon5OrderFaceCenter<weno5_JSchen>,
//     Recon5Order1DEulerEig<weno5_JSchen>, Recon5Order2DEulerEig<weno5_JSchen>,
//     Recon5OrderFaceCenter<weno5_Z>, Recon5Order1DEulerEig<weno5_Z>,
//     Recon5Order2DEulerEig<weno5_Z>, Recon5OrderFaceCenter<Teno5_Z>,
//     Recon5Order1DEulerEig<Teno5_Z>, Recon5Order2DEulerEig<Teno5_Z>,
//     Recon5OrderFaceCenter<Teno5_ZCT4>, Recon5Order1DEulerEig<Teno5_ZCT4>,
//     Recon5Order2DEulerEig<Teno5_ZCT4>, Recon5OrderFaceCenter<Teno5_ZCT7>,
//     Recon5Order1DEulerEig<Teno5_ZCT7>, Recon5Order2DEulerEig<Teno5_ZCT7>,
//     Recon5OrderFaceCenter<Teno5_CongZ>, Recon5Order1DEulerEig<Teno5_CongZ>,
//     Recon5Order2DEulerEig<Teno5_CongZ>,
//     Recon5OrderFaceCenter<Teno5_CongZCT4>,
//     Recon5Order1DEulerEig<Teno5_CongZCT4>,
//     Recon5Order2DEulerEig<Teno5_CongZCT4>,
//     Recon5OrderFaceCenter<Teno5_CongZCT7>,
//     Recon5Order1DEulerEig<Teno5_CongZCT7>,
//     Recon5Order2DEulerEig<Teno5_CongZCT7>>;

// using SolvePointFluxType =
//     std::variant<SolvePointSolver<3, fluxSolveEuler1D>,
//                  SolvePointSolver<1, fluxSolveLinearConv>,
//                  SolvePointSolver<1, fluxSolveBurgers>,
//                  SolvePointSolver<4, fluxSolveEuler2D>>;

// using FluxPointFluxType =
//     std::variant<FluxPointSolver<3, roeFlux1D2>, FluxPointSolver<3,
//     HLLCFlux1D>,
//                  FluxPointSolver<4, roeFlux2DSym>,
//                  FluxPointSolver<4, roeFlux2D>,
//                  FluxPointSolver<4, HLLCFlux2D2>>;
// using DifferenceType = std::variant<ExplicitDif6, MidNodeAndNodeDif6,
// Differ>;