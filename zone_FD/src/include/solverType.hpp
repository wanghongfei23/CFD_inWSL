#include <variant>

#include "differ.hpp"
#include "fluxPointFlux.hpp"
#include "fluxSchemes.hpp"
#include "proxy.h"
#include "reconstructor.hpp"
#include "reconstructor5order.hpp"
#include "solvePointFlux.hpp"

/*for ProxyDataManipulator*/
PRO_DEF_MEM_DISPATCH(Solve, solve);
PRO_DEF_MEM_DISPATCH(Init, init);
PRO_DEF_MEM_DISPATCH(SetConstNorm, setConstNorm);
PRO_DEF_MEM_DISPATCH(GetData, getData);
PRO_DEF_MEM_DISPATCH(Check, check);

/*for ProxySolvePointSolver additionally*/
PRO_DEF_MEM_DISPATCH(SetNlNr, setNlNr);

/*for ProxyFluxPointSolver*/

/*for ProxyDiffer*/
PRO_DEF_MEM_DISPATCH(SetConstantH, setConstantH);

struct ProxyDataManipulator
    : pro::facade_builder ::add_convention<Solve, void()>::
          add_convention<Check, void()>::
              add_convention<Init, void(DataReader, std::shared_ptr<OneDBnd>, std::shared_ptr<OneDBnd>)>::
                  add_convention<SetConstNorm, void(std::array<real, 3>)>::
                      add_convention<GetData, std::shared_ptr<Data>()>::
                          support_copy<pro::constraint_level::nontrivial>::build { };

struct ProxySolvePointSolver
    : pro::facade_builder ::add_convention<Solve, void()>::
          add_convention<Check, void()>::
              add_convention<Init, void(DataReader, std::shared_ptr<OneDBnd>, std::shared_ptr<OneDBnd>)>::
                  add_convention<SetConstNorm, void(std::array<real, 3>)>::
                      add_convention<SetNlNr, void(int, int)>::
                          add_convention<GetData, std::shared_ptr<Data>()>::
                              support_copy<pro::constraint_level::nontrivial>::build { };

struct ProxyDiffer
    : pro::facade_builder ::
          add_convention<Check, void()>::
              add_convention<Init, void(std::shared_ptr<Data>, std::shared_ptr<Data>, DataReader)>::
                  add_convention<SetConstantH, void(real)>::
                      add_convention<Solve, void()>::
                          support_copy<pro::constraint_level::nontrivial>::build { };

struct ProxyFluxPointSolver
    : pro::facade_builder ::
          add_convention<Check, void()>::
              add_convention<Solve, void()>::
                  add_convention<Init, void(std::shared_ptr<Data>, int, int)>::
                      add_convention<SetConstNorm, void(std::array<real, 3>)>::
                          add_convention<GetData, std::shared_ptr<Data>()>::
                              support_copy<pro::constraint_level::nontrivial>::build { };

class SolverType {
public:
    SolverType(Info* info_);
    SolverType() {};
    SolverType(const SolverType&);
    pro::proxy<ProxyDataManipulator> reconer;
    pro::proxy<ProxySolvePointSolver> solvePointSolver;
    pro::proxy<ProxyFluxPointSolver> fluxPointSolver;
    pro::proxy<ProxyDiffer> differ;

private:
    Info* info;
    void initReconer();
    void initSolvePointSolver();
    void initFluxPointSolver();
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