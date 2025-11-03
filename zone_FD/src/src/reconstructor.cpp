#include "reconstructor.hpp"
#include "eigenSystem.hpp"
#include "interScheme.hpp"

// std::shared_ptr<Data> DataManipulater::solvePointFluxes(int nl, int nr) {
//   assert(nl < bndL->getN());
//   assert(nr < bndR->getN());
//   fluxes = std::make_shared<Data>(n + nl + nr, nvar);
//   auto iterFlux = fluxes->begin();
//   for (int i = nl - 1; i >= 0; i--) {
//     auto iterVarbl = (*bndL)(i);
//     fluxSolve(iterVarbl, iterFlux, norm);
//     iterFlux += nvar;
//   }

//   for (int i = 0; i < n; i++) {
//     auto iterVar = varsReader(i);
//     fluxSolve(iterVar, iterFlux, norm);
//     iterFlux += nvar;
//   }
//   for (int i = 0; i < nr; i++) {
//     auto iterVarbr = (*bndR)(i);
//     fluxSolve(iterVarbr, iterFlux, norm);
//     iterFlux += nvar;
//   }
//   assert(iterFlux == fluxes->end());
// }

void Recon1Order::initVarsR() {
  assert(bndL->getN() >= 1);
  assert(bndR->getN() >= 1);
  int nPointR = n + 1;
  data = std::make_shared<Data>(n * 2, nvar);
}
void Recon1Order::leftBnd() {
  auto iterVar = varsReader(0);
  std::copy((*bndL)(0), (*bndL)(0) + nvar, iter);
  std::copy(iterVar, iterVar + nvar, iter + nvar);
  iter += 2 * nvar;
}
void Recon1Order::internal() {
  for (int i = 0; i < n - 2; i++) {
    auto iterVar = varsReader(i);
    std::copy(iterVar, iterVar + nvar, iter);
    std::copy(iterVar + nvar, iterVar + 2 * nvar, iter + nvar);
  }
  iter += 2 * nvar;
}
void Recon1Order::rightBnd() {
  auto iterVar = varsReader(n - 1);
  std::copy(iterVar, iterVar + nvar, iter);
  std::copy((*bndR)(0), (*bndR)(0) + nvar, iter + nvar);
  iter += 2 * nvar;
}