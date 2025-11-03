#include "bnds.hpp"

std::array<std::shared_ptr<OneDBnd>, 2> Bnds::getOneDBnd(int idim, int i, int j)
{
    std::array<std::shared_ptr<OneDBnd>, 2> res;
    int index;
    if (idim == 1) {
        index = (i + j * iMax[2]) * 2;
        res[0] = oneDBnds.at(index);
        res[1] = oneDBnds.at(index + 1);
    } else if (idim == 2) {
        index = (iMax[1] * iMax[2] + i + j * iMax[2]) * 2;
        res[0] = oneDBnds.at(index);
        res[1] = oneDBnds.at(index + 1);
    } else if (idim == 3) {
        index = (iMax[1] * iMax[2] + iMax[0] * iMax[2] + i + j * iMax[1]) * 2;
        res[0] = oneDBnds.at(index);
        res[1] = oneDBnds.at(index + 1);
    }
    return res;
}

void Bnds::update()
{
    for (auto ibnd : oneDBnds) {
        ibnd->update();
    }
}