#pragma once
#include "macro.hpp"
#include "cgnslib.h"
#include "block.hpp"
#include "data.hpp"
#include "info.hpp"


class CgnsIO
{
    public:
    void BlockCgnsOutput(Block* block,Info* info);
    void solCgnsOutput(Data* data,Info* info);
    void oneDsolOutput(Info* info);
    private:
};