#include "dataManipulater.hpp"
DataManipulater::DataManipulater(DataReader varsReader_,
                                 std::shared_ptr<OneDBnd> bndL_,
                                 std::shared_ptr<OneDBnd> bndR_)
    : varsReader(varsReader_), bndL(bndL_), bndR(bndR_) {
  n = varsReader.getN();
  nvar = varsReader.getNVar();
}
void DataManipulater::init(DataReader varsReader_,
                           std::shared_ptr<OneDBnd> bndL_,
                           std::shared_ptr<OneDBnd> bndR_) {
  varsReader = varsReader_;
  bndL = bndL_;
  bndR = bndR_;
  n = varsReader.getN();
  nvar = varsReader.getNVar();
}

void DataManipulater::setConstNorm(std::array<real, 3> norm_) { norm = norm_; }