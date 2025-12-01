#include "sp_distributor.hpp"

void SpDistributor::rhsSolve()
{
    if (dim <= 0 || dim > 3) {
        std::cout << "SpDistributor error: dim error";
        return;
    }
    if (dim >= 1) {
        long long timep = 0;
        SolverType solverTemp(solverType);
#pragma omp parallel for collapse(2) reduction(+ : timep) schedule(static) firstprivate(solverTemp)
        for (int i = 0; i < iMax[1]; i++)
            for (int j = 0; j < iMax[2]; j++) {
                constexpr int idim = 0;
                constexpr std::array<real, 3> normTemp = { 1.0, 0.0, 0.0 };
                auto oneDBnds = bnds->getOneDBnd(1, i, j);
                auto offsets = calOffset(1, i, j, iMax);

                DataReader primReader(iMax[idim], offsets[0], offsets[1], idim, prim);
                DataReader rhsReader(iMax[idim], offsets[0], offsets[1], idim, rhs);

                solverTemp.reconer->init(primReader, oneDBnds[0], oneDBnds[1]);
                solverTemp.reconer->setConstNorm(normTemp);
                solverTemp.reconer->solve();
                auto reconData=solverTemp.reconer->getData();
                int nData=reconData->getN();
                int nVarData=reconData->getN();

                solverTemp.fluxPointSolver->init(solverTemp.reconer->getData(), idim, prim->getNVar());
                solverTemp.fluxPointSolver->setConstNorm(normTemp);
                solverTemp.fluxPointSolver->solve();

                solverTemp.solvePointSolver->init(primReader, oneDBnds[0], oneDBnds[1]);
                solverTemp.solvePointSolver->setConstNorm(normTemp);
                solverTemp.solvePointSolver->solve();

                solverTemp.differ->init(solverTemp.fluxPointSolver->getData(), solverTemp.solvePointSolver->getData(), rhsReader);
                solverTemp.differ->setConstantH(info->geth(idim));
                solverTemp.differ->solve();
                // timep += spDis.timep;
            }
        timepp += timep;
    }

    if (dim >= 2) {
        long long timep = 0;
        SolverType solverTemp1(solverType);
#pragma omp parallel for collapse(2) reduction(+ : timep) schedule(static) firstprivate(solverTemp1)
        for (int i = 0; i < iMax[0]; i++)
            for (int j = 0; j < iMax[2]; j++) {
                constexpr int idim = 1;
                constexpr std::array<real, 3> normTemp = { 0.0, 1.0, 0.0 };

                auto oneDBnds = bnds->getOneDBnd(2, i, j);
                auto offsets = calOffset(2, i, j, iMax);

                DataReader primReader(iMax[idim], offsets[0], offsets[1], idim, prim);
                DataReader rhsReader(iMax[idim], offsets[0], offsets[1], idim, rhs);

                solverTemp1.reconer->init(primReader, oneDBnds[0], oneDBnds[1]);
                solverTemp1.reconer->setConstNorm(normTemp);
                solverTemp1.reconer->solve();

                solverTemp1.fluxPointSolver->init(solverTemp1.reconer->getData(), idim, prim->getNVar());
                solverTemp1.fluxPointSolver->setConstNorm(normTemp);
                solverTemp1.fluxPointSolver->solve();

                solverTemp1.solvePointSolver->init(primReader, oneDBnds[0], oneDBnds[1]);
                solverTemp1.solvePointSolver->setConstNorm(normTemp);
                solverTemp1.solvePointSolver->solve();

                solverTemp1.differ->init(solverTemp1.fluxPointSolver->getData(), solverTemp1.solvePointSolver->getData(), rhsReader);
                solverTemp1.differ->setConstantH(info->geth(idim));
                solverTemp1.differ->solve();
            }
        timepp += timep;
    }
}
