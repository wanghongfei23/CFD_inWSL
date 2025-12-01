// #pragma once
// #include "fluxSolver.hpp"
// #include "reconstructor.hpp"
// #include <concepts>

// template <typename T>
// concept RiemannSolver1DEuler = requires(T f, int a) {
//   { f(a, a, a, a, a, a) } -> std::same_as<std::array<real, 3>>;
// };

// template <typename T>
// concept FluxSolver1DEuler = requires(T f, int a) {
//   { f(a, a, a) } -> std::same_as<std::array<real, 3>>;
// };
// template <RiemannSolver1DEuler Solver>
// class FluxPointSolverEuler1D : public FluxPointSolver {
// public:
//   void solve() override;
// };

// template <FluxSolver1DEuler Flux>
// class SolvePointSolverEuler1D : public SolvePointSolver {
// public:
//   void initVarsR() override;
//   void leftBnd() override;
//   void internal() override;
//   void rightBnd() override;
// };