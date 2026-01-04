#include "sp_distributor.hpp"

void SpDistributor::rhsSolve()
{
    if (dim <= 0 || dim > 3) {
        std::cout << "SpDistributor error: dim error";
        return;
    }

    if (dim >= 1) {
        long long timep = 0;
#pragma omp parallel for collapse(2) reduction(+ : timep) schedule(static)
        for (int i = 0; i < iMax[1]; i++)
            for (int j = 0; j < iMax[2]; j++) {
                // std::cout<<std::format("thread num: {} \n",omp_get_thread_num());
                auto oneDBnds = bnds->getOneDBnd(1, i, j);
                auto offsets = calOffset(1, i, j, iMax);
                SpaceDis spDis(iMax[0], prim, rhs, oneDBnds[0], oneDBnds[1], info);

                spDis.setConstNorm({ 1, 0, 0 });
                spDis.setOffset(offsets[0], offsets[1]);
                spDis.setIDim(0);
                spDis.difference();
                timep += spDis.timep;
            }
        timepp += timep;
    }

    if (dim >= 2) {
        long long timep = 0;
#pragma omp parallel for collapse(2) reduction(+ : timep) schedule(static)
        for (int i = 0; i < iMax[0]; i++)
            for (int j = 0; j < iMax[2]; j++) {
                auto oneDBnds = bnds->getOneDBnd(2, i, j);
                auto offsets = calOffset(2, i, j, iMax);
                SpaceDis spDis(iMax[1], prim, rhs, oneDBnds[0], oneDBnds[1], info);

                spDis.setConstNorm({ 0, 1, 0 });
                spDis.setIDim(1);
                spDis.setOffset(offsets[0], offsets[1]);
                spDis.difference();
                timep += spDis.timep;
            }
        timepp += timep;
    }

    if (dim >= 3) {
        long long timep = 0;
#pragma omp parallel for collapse(2) reduction(+ : timep) schedule(static)
        for (int i = 0; i < iMax[0]; i++)
            for (int j = 0; j < iMax[1]; j++) {
                auto oneDBnds = bnds->getOneDBnd(3, i, j);
                auto offsets = calOffset(3, i, j, iMax);
                SpaceDis spDis(iMax[2], prim, rhs, oneDBnds[0], oneDBnds[1], info);

                spDis.setConstNorm({ 0, 0, 1 });
                spDis.setIDim(2);
                spDis.setOffset(offsets[0], offsets[1]);

                spDis.difference();
                timep += spDis.timep;
            }
        timepp += timep;
    }
}
