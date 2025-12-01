#pragma once
#include "dataManipulater.hpp"

class Reconstuctor : public DataManipulater {
public:
    Reconstuctor() {};
    Reconstuctor(DataReader varsR_, std::shared_ptr<OneDBnd> bndL_,
        std::shared_ptr<OneDBnd> bndR_)
        : DataManipulater(varsR_, bndL_, bndR_) {};
    void check() { std::cout << "initialized successfully Reconstuctor\n"; }
    void solve() final
    {
        initVarsR();
        iter = data->begin();
        leftBnd();
        internal();
        rightBnd();
        assert(iter == data->end());
    }
};

class Recon1Order : public Reconstuctor {
public:
protected:
    void initVarsR() override;
    void leftBnd() override;
    void internal() override;
    void rightBnd() override;
};
