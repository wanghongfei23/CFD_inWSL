// #include "EulerEquation1D.hpp"
// template <RiemannSolver1DEuler Solver>
// void FluxPointSolverEuler1D<Solver>::solve() {

//   assert(fluxes->getN() == valsR->getN() / 2);
//   auto iterFlux = fluxes->begin();
//   for (auto ivar = valsR->begin(); ivar != valsR->end(); ivar += 6) {
//     assert(ivar <= valsR->end());
//     // auto iflux = solver(ivar, norm);
//     // std::copy(iflux.begin(), iflux.end(), iterFlux);
//     iterFlux += 3;
//   }
// }