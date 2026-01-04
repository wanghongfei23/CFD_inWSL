

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
    WHFTCNSASF002,    // 18
    WHFTCNSAH002,    // 19
    WHFTCNSASF102,     // 20
    WHFTCNSASF103,     // 21
    WHFTCNSASF102_reciprocal,     // 22
    WHFTCNSASF103_reciprocal,     // 23
    WHFTCNSAS_fx,     // 24
    WHFTCNSAS_initial,     // 25
    WHFTCNSAS_approx_1,     // 26
    WHFTCNSAS_fx_real,     // 27
    WHFTCNSAS_approx_2,     // 28
    WHFTCNSASF202_2S,     // 29
    WHFTCNSASF202_NoS,     // 30
    WHFTCNSASF203_NoS,     // 31
    temp011,     // 32
    temp012,     // 33
    temp013,     // 34
    temp014,     // 35
    temp015,     // 36
    temp016,     // 37
    temp017,     // 38
    temp018,     // 39
    temp019     // 40
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