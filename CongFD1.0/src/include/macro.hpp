

#pragma once

#define real double
#define ind int

// #include <cmath>
// #include <stdio.h>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <map>

#include <chrono>
inline long timepp = 0;
inline long timesss = 0;

// 【王鸿飞】注意gamma
#define GAMMA 1.4   //其他
// #define GAMMA 5.0/3.0  //RTI

enum BndType {
    TYPENULL,
    PERIODIC1D,
    SYMMETRY1D,
    DIRICLET,
    DIRICLET_SODL,
    DIRICLET_SODR,
    FLUXGHOST,
    SUPERSONICOUTLET,
    SYMMETRYX, // only for 2D
    SYMMETRYY, // only for 2D
    DoubleMachUp,
    ARDBND
};
enum InterMethod {
    FIRSTORDER,      // 0
    MUSCL,           // 1
    WCNS5,           // 2
    WCNSZ5,          // 3
    WCNS5Char,       // 4
    WCNSZ5Char,      // 5
    WCNS5CONG,       // 6
    TCNSCongA,       // 7
    WCNS5CONGZ,      // 8
    WCNS5CONGZCT4,   // 9
    WCNS5CONGZCT7,   // 10
    TCNS5,           // 11
    TCNS5CT4,        // 12
    TCNS5CT7,        // 13
    LINEAR5,         // 14
    NICEST5,         // 15
    INTERMAX,        // 16
    WHFTCNSA,        // 17
    WHFTCNSAF002,    // 18
    WHFTCNSAH002,    // 19
    WHFTCNSAF102,     // 20
    ending     // 21
};

enum DiffMethod {
    HDS6,
    TRAD6,
    TRAD2,
    MND6
};

enum EquationType {
    LINEARCONV1D,
    BURGERS1D,
    EULER,
    ACCURACYTEST
};

enum FluxMethod {
    HLLC1D,
    ROE1D,
    HLLC2D
};

enum TimeMethod {
    RK3SSP,
    EulerFront
};

enum SourceType {
    SOURCENULL,
    GRAVITY
};

// int index(int,int,int,std::array<int,3>);
// std::array<int,2> calOffset(int dim,int i,int j,std::array<int,3>);
// std::array<int,2> calOffsetInverse(int idim,int i,int j,std::array<int,3> iMax);

constexpr int index(int i, int j, int k, std::array<int, 3> iMax)
{
    return i + j * iMax[0] + k * iMax[0] * iMax[1];
}

constexpr std::array<int, 2> calOffset(int idim, int i, int j, std::array<int, 3> iMax)
{
    std::array<int, 3> offsets { 1, iMax[0], iMax[0] * iMax[1] };
    std::array<int, 2> res;
    if (idim == 1) {
        res[0] = i * offsets[1] + j * offsets[2]; // i0
        res[1] = offsets[0]; // offset
    } else if (idim == 2) {
        res[0] = i * offsets[0] + j * offsets[2]; // i0
        res[1] = offsets[1]; // offset
    } else if (idim == 3) {
        res[0] = i * offsets[0] + j * offsets[1]; // i0
        res[1] = offsets[2]; // offset
    }
    return res;
}

constexpr std::array<int, 2> calOffsetInverse(int idim, int i, int j, std::array<int, 3> iMax)
{
    std::array<int, 3> offsets { 1, iMax[0], iMax[0] * iMax[1] };
    std::array<int, 2> res;
    if (idim == 1) {
        res[0] = i * offsets[1] + j * offsets[2] + (iMax[0] - 1) * offsets[0]; // i0
        res[1] = -offsets[0]; // offset
    } else if (idim == 2) {
        res[0] = i * offsets[0] + j * offsets[2] + (iMax[1] - 1) * offsets[1]; // i0
        res[1] = -offsets[1]; // offset
    } else if (idim == 3) {
        res[0] = i * offsets[0] + j * offsets[1] + (iMax[2] - 1) * offsets[2]; // i0
        res[1] = -offsets[2]; // offset
    }
    return res;
}