#include "solverType.hpp"
#include <variant>

int main() {
  auto bndl = std::make_shared<OneDBnd>(10, 3, SUPERSONICOUTLET);
  auto bndr = std::make_shared<OneDBnd>(10, 3, SUPERSONICOUTLET);
  auto data = new Data(50, 3);
  std::vector<real> initvalues(150);
  for (int i = 0; i < 150; i++) {
    initvalues[i] = i;
  }
  data->setValue(initvalues);

  DataReader dataReader(50, 0, 1, 0, data);
  // = std::make_shared<DataReader>(50, 0, 1, 0, data);
  Info *info = new Info;

  // using solvePointFluxSolvers =
  //     std::variant<SolvePointSolver<3, fluxSolveEuler1D>,
  //                  SolvePointSolver<4, fluxSolveEuler2D>>;
  // auto func = [=](int ii, DataReader varsReader_,
  //                 std::shared_ptr<OneDBnd> bndL_,
  //                 std::shared_ptr<OneDBnd> bndR_, Info *info_, int nl_,
  //                 int nr_) -> solvePointFluxSolvers {
  //   if (ii == 3)
  //     return SolvePointSolver<3, fluxSolveEuler1D>(dataReader, bndl, bndr,
  //     info,
  //                                                  1, 1);
  //   else
  //     return SolvePointSolver<4, fluxSolveEuler2D>(dataReader, bndl, bndr,
  //     info,
  //                                                  1, 1);
  // };

  Recon5OrderFaceCenter<weno5_JSchen> Reconer(dataReader, bndl, bndr);
  SolvePointSolver<3, fluxSolveEuler1D> solvePointFlux(dataReader, bndl, bndr,
                                                       1, 1);
  solvePointFlux.setConstNorm({1, 0, 0});

  Reconer.solve();
  FluxPointSolver<3, roeFlux1D2> fluxPointFlux(Reconer.data, 0, 3);
  fluxPointFlux.setConstNorm({1, 0, 0});
  solvePointFlux.solve();
  fluxPointFlux.solve();
  SolverType aa(info);
  aa.reconer->init(dataReader, bndl, bndr);

  return 0;
}