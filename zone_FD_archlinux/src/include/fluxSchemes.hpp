#pragma once

#include "macro.hpp"
#include <algorithm>
#include <array>
#include <span>
#include <vector>

// 一维Roe通量计算相关类型定义
typedef std::array<real, 2> arr2;

// 一维Roe通量计算函数（第二种实现）
inline void roeFlux1D2(const std::span<real, 6>& iter,
    const std::span<real, 3>& res,
    const std::array<real, 3>& norm)
{
    // reference: https://blog.csdn.net/Tankrun1997/article/details/132743487
    real rl = iter[0];
    real ul = iter[1];
    real pl = iter[2];
    real rr = iter[3];
    real ur = iter[4];
    real pr = iter[5];
    enum { L,
        R };

    real gamma = GAMMA;

    arr2 H = { pl / rl * GAMMA / (GAMMA - 1) + (ul * ul) / 2,
        pr / rr * GAMMA / (GAMMA - 1) + (ur * ur) / 2 };

    real rBar = std::sqrt(rl * rr);
    real uBar = (ul * std::sqrt(rl) + ur * std::sqrt(rr)) / (std::sqrt(rl) + std::sqrt(rr));
    real HBar = (H[L] * std::sqrt(rl) + H[R] * std::sqrt(rr)) / (std::sqrt(rl) + std::sqrt(rr));
    real cBar = std::sqrt((gamma - 1) * (HBar - uBar * uBar / 2));
    real cBar2 = (gamma - 1) * (HBar - uBar * uBar / 2);

    real dr = rr - rl;
    real du = ur - ul;
    real dp = pr - pl;

    real K1[3] = { 1, uBar, uBar * uBar / 2 };
    real K2[3] = { 1, uBar - cBar, HBar - uBar * cBar };
    real K3[3] = { 1, uBar + cBar, HBar + uBar * cBar };

    real lambda[3] = { std::abs(uBar), std::abs(uBar + cBar),
        std::abs(uBar - cBar) };

    real eps = 0.1 * (std::abs(uBar) + cBar);
    for (int i = 0; i < 3; i++) {
        if (lambda[i] < eps) {
            lambda[i] = (lambda[i] * lambda[i] + eps * eps) / (2.0 * eps);
        }
    }

    real alpha[5] = { lambda[0] * (dr - dp / cBar2),
        lambda[1] / (2.0 * cBar2) * (dp + rBar * cBar * du),
        lambda[2] / (2.0 * cBar2) * (dp - rBar * cBar * du) };
    alpha[3] = alpha[0] + alpha[1] + alpha[2];
    alpha[4] = cBar * (alpha[1] - alpha[2]);

    real FL[3] = { rl * ul, rl * ul * ul + pl, rl * H[L] * ul };
    real FR[3] = { rr * ur, rr * ur * ur + pr, rr * H[R] * ur };

    for (int i = 0; i < 3; i++) {
        res[i] = (FL[i] + FR[i]) / 2.0;
    }
    res[0] -= 0.5 * alpha[3];
    res[1] -= 0.5 * (alpha[3] * uBar + alpha[4]);
    res[2] -= 0.5 * (HBar * alpha[3] + uBar * alpha[4] - cBar2 * alpha[0] / (gamma - 1));
}

// 一维HLLC通量计算函数
inline void HLLCFlux1D(const std::span<real, 6>& iter,
    const std::span<real, 3>& res,
    const std::array<real, 3>& norm)
{
    // reference:https://zhuanlan.zhihu.com/p/583555029
    real rl = iter[0];
    real ul = iter[1];
    real pl = iter[2];
    real rr = iter[3];
    real ur = iter[4];
    real pr = iter[5];
    enum { L,
        R };
    arr2 H = { pl / rl * GAMMA / (GAMMA - 1) + (ul * ul) / 2,
        pr / rr * GAMMA / (GAMMA - 1) + (ur * ur) / 2 };
    arr2 RT = { pl / rl, pr / rr };
    real gamma = GAMMA;
    real cl = std::sqrt(gamma * RT[L]);
    real cr = std::sqrt(gamma * RT[R]);

    real uBar = (ul + ur) / 2, cBar = (cl + cr) / 2;
    real SL = std::min(ur - cr, ul - cl);
    real SR = std::max(ul + cl, ur + cr);

    real Sstar = ((pr - pl) + (rl * ul * (SL - ul) - rr * ur * (SR - ur))) / (rl * (SL - ul) - rr * (SR - ur));
    if (SL >= 0) {
        res[0] = rl * ul;
        res[1] = rl * ul * ul + pl;
        res[2] = rl * H[L] * ul;
    } else if (Sstar >= 0) {
        real pStar = pl + rl * (SL - ul) * (Sstar - ul);
        real U[3] = { rl, rl * ul, pl / (gamma - 1) + rl * ul * ul / 2 };
        real F[3] = { rl * ul, rl * ul * ul + pl, rl * H[L] * ul };
        real D[3] = { 0, 1, Sstar };
        for (int i = 0; i < 3; i++) {
            res[i] = (Sstar * (SL * U[i] - F[i]) + SL * pStar * D[i]) / (SL - Sstar);
        }

    } else if (SR >= 0) {
        real pStar = pr + rr * (SR - ur) * (Sstar - ur);
        real U[3] = { rr, rr * ur, pr / (gamma - 1) + rr * ur * ur / 2 };
        real F[3] = { rr * ur, rr * ur * ur + pr, rr * H[R] * ur };
        real D[3] = { 0, 1, Sstar };
        for (int i = 0; i < 3; i++) {
            res[i] = (Sstar * (SR * U[i] - F[i]) + SR * pStar * D[i]) / (SR - Sstar);
        }
    } else {
        res[0] = rr * ur;
        res[1] = rr * ur * ur + pr;
        res[2] = rr * H[R] * ur;
    }
}

// 二维HLLC通量计算函数
inline void HLLCFlux2D(const std::span<real, 8>& iter,
    const std::span<real, 4>& res,
    const std::array<real, 3>& norm)
{

    real rl = iter[0];
    real ul = iter[1];
    real vl = iter[2];
    real pl = iter[3];
    real rr = iter[4];
    real ur = iter[5];
    real vr = iter[6];
    real pr = iter[7];
    enum { L,
        R };
    // reference:https://zhuanlan.zhihu.com/p/583555029
    // if(pl<0 || pr<0 || rr<0 || rl<0)
    // {
    //     std::cout<<"fluxScheme error: positivity break\n";
    // }

    arr2 H;
    H[L] = (ul * ul + vl * vl) / 2 + pl / rl * GAMMA / (GAMMA - 1);
    H[R] = (ur * ur + vr * vr) / 2 + pr / rr * GAMMA / (GAMMA - 1);

    arr2 Vn = { (norm[1] > norm[0]) ? vl : ul, (norm[1] > norm[0]) ? vr : ur };
    real gamma = GAMMA;
    real cl = std::sqrt((gamma - 1) * pl / rl * GAMMA / (GAMMA - 1));
    real cr = std::sqrt((gamma - 1) * pr / rr * GAMMA / (GAMMA - 1));

    real uBar = (Vn[L] + Vn[R]) / 2, cBar = (cl + cr) / 2;
    // real SL=std::min(uBar-cBar,Vn[L]-cl);
    // real SR=std::max(uBar+cBar,Vn[R]+cr);

    real SL = std::min(Vn[R] - cr, Vn[L] - cl);
    real SR = std::max(Vn[L] + cl, Vn[R] + cr);

    real Sstar = ((pr - pl) + (rl * Vn[L] * (SL - Vn[L]) - rr * Vn[R] * (SR - Vn[R]))) / (rl * (SL - Vn[L]) - rr * (SR - Vn[R]));
    if (SL > 0) {
        res[0] = rl * Vn[L];
        res[1] = ul * Vn[L] * rl + pl * norm[0];
        res[2] = vl * Vn[L] * rl + pl * norm[1];
        res[3] = rl * H[L] * Vn[L];
    } else if (Sstar > 0) {
        real pStar = pl + rl * (SL - Vn[L]) * (Sstar - Vn[L]);
        real U[4] = { rl, rl * ul, rl * vl, rl * H[L] - pl };
        real F[4] = { rl * Vn[L], ul * Vn[L] * rl + pl * norm[0],
            vl * Vn[L] * rl + pl * norm[1], rl * H[L] * Vn[L] };
        real D[4] = { 0, norm[0], norm[1], Sstar };
        for (int i = 0; i < 4; i++) {
            res[i] = (SL * (Sstar * U[i] + pStar * D[i]) - Sstar * F[i]) / (SL - Sstar);
        }

    } else if (SR > 0) {
        real pStar = pr + rr * (SR - Vn[R]) * (Sstar - Vn[R]);
        real U[4] = { rr, rr * ur, rr * vr, rr * H[R] - pr };

        real F[4] = { rr * Vn[R], ur * Vn[R] * rr + pr * norm[0],
            vr * Vn[R] * rr + pr * norm[1], rr * H[R] * Vn[R] };

        real D[4] = { 0, norm[0], norm[1], Sstar };
        for (int i = 0; i < 4; i++) {
            res[i] = (SR * (Sstar * U[i] + pStar * D[i]) - Sstar * F[i]) / (SR - Sstar);
        }

    } else {
        res[0] = rr * Vn[R];
        res[1] = ur * Vn[R] * rr + pr * norm[0];
        res[2] = vr * Vn[R] * rr + pr * norm[1];
        res[3] = rr * H[R] * Vn[R];
    }
}

// 二维HLLC通量计算函数（第二种实现）
inline void HLLCFlux2D2(const std::span<real, 8>& iter,
    const std::span<real, 4>& res,
    const std::array<real, 3>& norm)
{
    //chenyuqing: HLLC的对称形式

    real rl = iter[0];
    real ul = iter[1];
    real vl = iter[2];
    real pl = iter[3];
    real rr = iter[4];
    real ur = iter[5];
    real vr = iter[6];
    real pr = iter[7];
    enum { L,
        R };
    // if(pl<0 || pr<0 || rr<0 || rl<0)
    // {
    //     std::cout<<"fluxScheme error: positivity break\n";
    // }

    arr2 H;
    H[L] = (ul * ul + vl * vl) / 2 + pl / rl * GAMMA / (GAMMA - 1);
    H[R] = (ur * ur + vr * vr) / 2 + pr / rr * GAMMA / (GAMMA - 1);
    //chenyuqing：修改HLLC焓计算顺序
    // H[L] = pl / rl * GAMMA / (GAMMA - 1) + ul * ul /2 + vl * vl / 2 ;
    // H[R] = pr / rr * GAMMA / (GAMMA - 1) + ur * ur /2 + vr * vr / 2 ;
    arr2 Vn = { (norm[1] > norm[0]) ? vl : ul, (norm[1] > norm[0]) ? vr : ur };
    real gamma = GAMMA;
    real cl = std::sqrt(pl / rl * gamma);
    real cr = std::sqrt(pr / rr * gamma);

    real rBar = (rl + rr) / 2, cBar = (cl + cr) / 2;
    real ps = std::max(0.0, (pl + pr) / 2 - rBar * cBar * (Vn[R] - Vn[L]) / 2);
    real ql = (ps <= pl)
        ? 1.0
        : std::sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (ps / pl - 1.0));
    real qr = (ps <= pr)
        ? 1.0
        : std::sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (ps / pr - 1.0));
    real SL = Vn[L] - cl * ql, SR = Vn[R] + cr * qr;

    // real uBar=(Vn[L]+Vn[R])/2,cBar=(cl+cr)/2;
    // real SL=std::min(uBar-cBar,Vn[L]-cl);
    // real SR=std::max(uBar+cBar,Vn[R]+cr);

    // real SL=std::min(Vn[R]-cr,Vn[L]-cl);
    // real SR=std::max(Vn[L]+cl,Vn[R]+cr);

    real Sstar = ((pr - pl) + (rl * Vn[L] * (SL - Vn[L]) - rr * Vn[R] * (SR - Vn[R]))) / (rl * (SL - Vn[L]) - rr * (SR - Vn[R]));
    //chenyuqing：修改HLLC中间波速公式
    //real Sstar = (pr - pl + rl * Vn[L] * (SL - Vn[L]) - rr * Vn[R] * (SR - Vn[R])) / (rl * SL + rr * Vn[R] - rl * Vn[L] -  rr * SR );
    if (SL > 0) {
        res[0] = rl * Vn[L];
        res[1] = rl * ul * Vn[L] + pl * norm[0];
        res[2] = rl * vl * Vn[L] + pl * norm[1];
        res[3] = rl * H[L] * Vn[L];
    } else if (Sstar > 0) {
        real U[4] = { rl, rl * ul, rl * vl, rl * H[L] - pl };
        real F[4] = { rl * Vn[L], rl * ul * Vn[L] + pl * norm[0],
            rl * vl * Vn[L] + pl * norm[1], rl * H[L] * Vn[L] };
        real Usf = (SL - Vn[L]) / (SL - Sstar);
        real Ustar[4] = {
            Usf * rl, Usf * (Sstar * norm[0] + ul * norm[1]) * rl,
            Usf * (Sstar * norm[1] + vl * norm[0]) * rl,
            Usf * (U[3] + (Sstar - Vn[L]) * (rl * Sstar + pl / (SL - Vn[L])))
        };
        for (int i = 0; i < 4; i++) {
            res[i] = F[i] + SL * (Ustar[i] - U[i]);
        }
    } 
    //chenyuqing：修改HLLC通量函数取值
    else if (Sstar == 0) {
        {
            real U[4] = { rl, rl * ul, rl * vl, rl * H[L] - pl };
            real F[4] = { rl * Vn[L], rl * ul * Vn[L] + pl * norm[0],
                rl * vl * Vn[L] + pl * norm[1], rl * H[L] * Vn[L] };
            real Usf = (SL - Vn[L]) / (SL - Sstar);
            real Ustar[4] = {
                Usf * rl, Usf * (Sstar * norm[0] + ul * norm[1]) * rl,
                Usf * (Sstar * norm[1] + vl * norm[0]) * rl,
                Usf * (U[3] + (Sstar - Vn[L]) * (rl * Sstar + pl / (SL - Vn[L])))
            };
            for (int i = 0; i < 4; i++) {
                res[i] = (F[i] + SL * (Ustar[i] - U[i])) / 2;
            }
        }

        {
            real U[4] = { rr, rr * ur, rr * vr, rr * H[R] - pr };

            real F[4] = { rr * Vn[R], rr * ur * Vn[R] + pr * norm[0],
                rr * vr * Vn[R] + pr * norm[1], rr * H[R] * Vn[R] };
            real Usf = (SR - Vn[R]) / (SR - Sstar);
            real Ustar[4] = {
                Usf * rr, Usf * (Sstar * norm[0] + ur * norm[1]) * rr,
                Usf * (Sstar * norm[1] + vr * norm[0]) * rr,
                Usf * (U[3] + (Sstar - Vn[R]) * (rr * Sstar + pr / (SR - Vn[R])))
            };
            for (int i = 0; i < 4; i++) {
                res[i] += (F[i] + SR * (Ustar[i] - U[i])) / 2;
            }
        }

    }
     else if (SR > 0) {
        real U[4] = { rr, rr * ur, rr * vr, rr * H[R] - pr };

        real F[4] = { rr * Vn[R], rr * ur * Vn[R] + pr * norm[0],
            rr * vr * Vn[R] + pr * norm[1], rr * H[R] * Vn[R] };
        real Usf = (SR - Vn[R]) / (SR - Sstar);
        real Ustar[4] = {
            Usf * rr, Usf * (Sstar * norm[0] + ur * norm[1]) * rr,
            Usf * (Sstar * norm[1] + vr * norm[0]) * rr,
            Usf * (U[3] + (Sstar - Vn[R]) * (rr * Sstar + pr / (SR - Vn[R])))
        };
        for (int i = 0; i < 4; i++) {
            res[i] = F[i] + SR * (Ustar[i] - U[i]);
        }
    } else {
        res[0] = rr * Vn[R];
        res[1] = rr * ur * Vn[R] + pr * norm[0];
        res[2] = rr * vr * Vn[R] + pr * norm[1];
        res[3] = rr * H[R] * Vn[R];
    }
}

// 二维Roe通量计算函数
inline void roeFlux2D(const std::span<real, 8>& iter,
    const std::span<real, 4>& res,
    const std::array<real, 3>& norm)
{

    real rl = iter[0];
    real ul = iter[1];
    real vl = iter[2];
    real pl = iter[3];
    real rr = iter[4];
    real ur = iter[5];
    real vr = iter[6];
    real pr = iter[7];
    enum { L,
        R };
    std::array<real, 4> FcL, FcR;
    if (pl < 0 || pr < 0) {
        std::cout << "fluxScheme error: Pressure positivity break\n";
    }
    if (rr < 0 || rl < 0) {
        std::cout << "fluxScheme error: Density positivity break\n";
    }

    arr2 H;
    H[L] = (ul * ul + vl * vl) / 2 + pl / rl * GAMMA / (GAMMA - 1);
    H[R] = (ur * ur + vr * vr) / 2 + pr / rr * GAMMA / (GAMMA - 1);

    arr2 Vn = { ul * norm[0] + vl * norm[1], ur * norm[0] + vr * norm[1] };
    FcL[0] = rl * Vn[L];
    FcL[1] = ul * Vn[L] * rl + pl * norm[0];
    FcL[2] = vl * Vn[L] * rl + pl * norm[1];
    FcL[3] = rl * H[L] * Vn[L];

    FcR[0] = rr * Vn[R];
    FcR[1] = ur * Vn[R] * rr + pr * norm[0];
    FcR[2] = vr * Vn[R] * rr + pr * norm[1];
    FcR[3] = rr * H[R] * Vn[R];

    double rhoAvg, uAvg, vAvg, HAvg, cAvg, VnAvg, q_2Avg, coef1, coef2;
    coef1 = std::sqrt(rl);
    coef2 = std::sqrt(rr);
    real divisor = 1.0 / (std::sqrt(rl) + std::sqrt(rr));
    rhoAvg = std::sqrt(rl * rr);
    uAvg = (coef1 * ul + coef2 * ur) * divisor;
    vAvg = (coef1 * vl + coef2 * vr) * divisor;
    HAvg = (coef1 * H[L] + coef2 * H[R]) * divisor;
    q_2Avg = (uAvg * uAvg + vAvg * vAvg) / 2;

    cAvg = std::sqrt((GAMMA - 1) * (HAvg - q_2Avg));
    if ((HAvg - q_2Avg) < 0) {
        cAvg = std::sqrt(0.1);
    }
    VnAvg = uAvg * norm[0] + vAvg * norm[1];

    double lambda[3] = { std::abs(VnAvg - cAvg), std::abs(VnAvg),
        std::abs(VnAvg + cAvg) };
    real eps = 0.05 * (std::abs(VnAvg) + cAvg);
    for (int i = 0; i < 3; i++) {
        if (lambda[i] < eps) {
            lambda[i] = (lambda[i] * lambda[i] + eps * eps) / (2.0 * eps);
        }
    }

    double deltaP = pr - pl, deltaVn = Vn[R] - Vn[L], deltaU = ur - ul,
           deltaV = vr - vl, deltaRho = rr - rl, coef;
    double FDispassion[4];
    coef1 = (deltaP - rhoAvg * cAvg * deltaVn) / (2 * cAvg * cAvg);
    FDispassion[0] = lambda[0] * coef1 * 1;
    FDispassion[1] = lambda[0] * coef1 * (uAvg - cAvg * norm[0]);
    FDispassion[2] = lambda[0] * coef1 * (vAvg - cAvg * norm[1]);
    FDispassion[3] = lambda[0] * coef1 * (HAvg - cAvg * VnAvg);

    coef2 = deltaRho - deltaP / (cAvg * cAvg);
    FDispassion[0] += lambda[1] * (coef2 * 1.0 + rhoAvg * 0.0);
    FDispassion[1] += lambda[1] * (coef2 * uAvg + rhoAvg * (deltaU - deltaVn * norm[0]));
    FDispassion[2] += lambda[1] * (coef2 * vAvg + rhoAvg * (deltaV - deltaVn * norm[1]));
    FDispassion[3] += lambda[1] * (coef2 * q_2Avg + rhoAvg * (uAvg * deltaU + vAvg * deltaV - VnAvg * deltaVn));

    coef = (deltaP + rhoAvg * cAvg * deltaVn) / (2.0 * cAvg * cAvg);
    FDispassion[0] += lambda[2] * (coef * 1);
    FDispassion[1] += lambda[2] * (coef * (uAvg + cAvg * norm[0]));
    FDispassion[2] += lambda[2] * (coef * (vAvg + cAvg * norm[1]));
    FDispassion[3] += lambda[2] * (coef * (HAvg + cAvg * VnAvg));

    for (int i = 0; i < 4; i++) {
        res[i] = 0.5 * (FcL[i] + FcR[i] - FDispassion[i]);
    }
}

// 二维Roe通量计算函数（对称版本）
inline void roeFlux2DSym(const std::span<real, 8>& iter,
    const std::span<real, 4>& res,
    const std::array<real, 3>& norm)
{
    //chenyuqing: ROE的对称格式

    real rl = iter[0];
    real ul = iter[1];
    real vl = iter[2];
    real pl = iter[3];
    real rr = iter[4];
    real ur = iter[5];
    real vr = iter[6];
    real pr = iter[7];
    enum { L,
        R };
    std::array<real, 4> FcL, FcR;
    if (pl < 0 || pr < 0) {
        std::cout << "fluxScheme error: Pressure positivity break\n";
    }
    if (rr < 0 || rl < 0) {
        std::cout << "fluxScheme error: Density positivity break\n";
    }

    arr2 H;
    H[L] = (ul * ul + vl * vl) / 2 + pl / rl * GAMMA / (GAMMA - 1);
    H[R] = (ur * ur + vr * vr) / 2 + pr / rr * GAMMA / (GAMMA - 1);

    arr2 Vn = { (norm[0] > norm[1]) ? ul : vl, (norm[0] > norm[1]) ? ur : vr };
    FcL[0] = rl * Vn[L];
    FcL[1] = rl * ul * Vn[L] + pl * norm[0];
    FcL[2] = rl * vl * Vn[L] + pl * norm[1];
    //chenyuqing：修改Roe格式中左值通量函数表达式3
    //FcL[2] = rl * Vn[L] * vl + pl * norm[1];
    FcL[3] = rl * H[L] * Vn[L];

    FcR[0] = rr * Vn[R];
    FcR[1] = rr * ur * Vn[R] + pr * norm[0];
    FcR[2] = rr * vr * Vn[R] + pr * norm[1];
    //chenyuqing：修改Roe格式中右值通量函数表达式3
    //FcR[2] = rr * Vn[R] * vr + pr * norm[1];
    FcR[3] = rr * H[R] * Vn[R];

    double rhoAvg, uAvg, vAvg, HAvg, cAvg, VnAvg, q_2Avg, coef1, coef2;
    coef1 = std::sqrt(rl);
    coef2 = std::sqrt(rr);
    real divisor = 1.0 / (std::sqrt(rl) + std::sqrt(rr));

    rhoAvg = std::sqrt(rl * rr);
    uAvg = (coef1 * ul + coef2 * ur) * divisor;
    vAvg = (coef1 * vl + coef2 * vr) * divisor;
    HAvg = (coef1 * H[L] + coef2 * H[R]) * divisor;
    q_2Avg = (uAvg * uAvg + vAvg * vAvg) / 2;
    cAvg = std::sqrt((GAMMA - 1) * (HAvg - q_2Avg));

    //chenyuqing：修改Roe平均值
    // real co1 = std::sqrt(rl) / (std::sqrt(rl) + std::sqrt(rr));
    // real co2 = 1.0 - co1;
    // rhoAvg = std::sqrt(rl * rr);
    // uAvg = co1 * ul + co2 * ur;
    // vAvg = co1 * vl + co2 * vr;
    // HAvg = co1 * H[L] + co2 * H[R];
    // q_2Avg = (uAvg * uAvg + vAvg * vAvg) / 2;
    // cAvg = std::sqrt((GAMMA - 1) * (HAvg - q_2Avg));

    // if((HAvg-q_2Avg)<0)
    // {
    //     cAvg=std::sqrt(0.1);
    // }
    VnAvg = (norm[0] > norm[1]) ? uAvg : vAvg;

    double lambda[3] = { std::abs(VnAvg - cAvg), std::abs(VnAvg),
        std::abs(VnAvg + cAvg) };
    // real eps=0.1*(std::abs(VnAvg)+cAvg);
    // for(int i=0;i<3;i++)
    // {
    //     if(lambda[i]<eps)
    //     {
    //         lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
    //     }
    // }
    real eps = 0.05 * (std::abs(VnAvg) + cAvg);
    for (int i = 0; i < 3; i++) {
        if (lambda[i] < eps) {
            lambda[i] = (lambda[i] * lambda[i] + eps * eps) / (2.0 * eps);
        }
    }

    double deltaP = pr - pl, deltaVn = Vn[R] - Vn[L], deltaU = ur - ul,
           deltaV = vr - vl, deltaRho = rr - rl, coef;
    double FDispassion[4];
    //chenyuqing:正确求和顺序
    coef1 = (deltaP - rhoAvg * cAvg * deltaVn) / (2.0 * cAvg * cAvg);
    FDispassion[0] = lambda[0] * (coef1 * 1);
    FDispassion[1] = lambda[0] * (coef1 * (uAvg - cAvg * norm[0]));
    FDispassion[2] = lambda[0] * (coef1 * (vAvg - cAvg * norm[1]));
    FDispassion[3] = lambda[0] * (coef1 * (HAvg - cAvg * VnAvg));

    coef = (deltaP + rhoAvg * cAvg * deltaVn) / (2.0 * cAvg * cAvg);
    FDispassion[0] += lambda[2] * (coef * 1);
    FDispassion[1] += lambda[2] * (coef * (uAvg + cAvg * norm[0]));
    FDispassion[2] += lambda[2] * (coef * (vAvg + cAvg * norm[1]));
    FDispassion[3] += lambda[2] * (coef * (HAvg + cAvg * VnAvg));

    coef2 = deltaRho - deltaP / (cAvg * cAvg);
    FDispassion[0] += lambda[1] * (coef2 * 1.0 + rhoAvg * 0.0);
    FDispassion[1] += lambda[1] * (coef2 * uAvg + rhoAvg * ((norm[1] > norm[0]) ? deltaU : 0));
    FDispassion[2] += lambda[1] * (coef2 * vAvg + rhoAvg * ((norm[1] > norm[0]) ? 0 : deltaV));
    FDispassion[3] += lambda[1] * (coef2 * q_2Avg + rhoAvg * ((norm[1] > norm[0]) ? uAvg * deltaU : vAvg * deltaV));


    //chenyuqing：修改Roe特征值求和顺序
    // coef1 = (deltaP - rhoAvg * cAvg * deltaVn) / (2.0 * cAvg * cAvg);
    // FDispassion[0] = lambda[0] * (coef1 * 1.0);
    // FDispassion[1] = lambda[0] * (coef1 * (uAvg - cAvg * norm[0]));
    // FDispassion[2] = lambda[0] * (coef1 * (vAvg - cAvg * norm[1]));
    // FDispassion[3] = lambda[0] * (coef1 * (HAvg - cAvg * VnAvg));

    // coef2 = deltaRho - deltaP / (cAvg * cAvg);
    // FDispassion[0] += lambda[1] * (coef2 * 1.0 + rhoAvg * 0.0);
    // FDispassion[1] += lambda[1] * (coef2 * uAvg + rhoAvg * ((norm[1] > norm[0]) ? deltaU : 0));
    // FDispassion[2] += lambda[1] * (coef2 * vAvg + rhoAvg * ((norm[1] > norm[0]) ? 0 : deltaV));
    // FDispassion[3] += lambda[1] * (coef2 * q_2Avg + rhoAvg * ((norm[1] > norm[0]) ? uAvg * deltaU : vAvg * deltaV));

    // coef = (deltaP + rhoAvg * cAvg * deltaVn) / (2.0 * cAvg * cAvg);
    // FDispassion[0] += lambda[2] * (coef * 1.0);
    // FDispassion[1] += lambda[2] * (coef * (uAvg + cAvg * norm[0]));
    // FDispassion[2] += lambda[2] * (coef * (vAvg + cAvg * norm[1]));
    // FDispassion[3] += lambda[2] * (coef * (HAvg + cAvg * VnAvg));

    for (int i = 0; i < 4; i++) {
        res[i] = (FcL[i] + FcR[i] - FDispassion[i]) / 2;
    }
}

// 线性对流方程通量计算函数
inline void fluxSolveLinearConv(const std::span<real, 1>& iterVar,
    const std::span<real, 1>& iterFlux,
    const std::array<real, 3>& norm)
{
    iterFlux[0] = iterVar[0];
}

// Burgers方程通量计算函数
inline void fluxSolveBurgers(const std::span<real, 1>& iterVar,
    const std::span<real, 1>& iterFlux,
    const std::array<real, 3>& norm)
{
    iterFlux[0] = iterVar[0] * iterVar[0];
}

// 一维欧拉方程通量计算函数
inline void fluxSolveEuler1D(const std::span<real, 3>& iterVar,
    const std::span<real, 3>& iterFlux,
    const std::array<real, 3>& norm)
{
    real r = iterVar[0], u = iterVar[1], p = iterVar[2];
    iterFlux[0] = r * u;
    iterFlux[1] = r * u * u + p;
    iterFlux[2] = (GAMMA / (GAMMA - 1) * p + r * (u * u) / 2) * u;
}

// 二维欧拉方程通量计算函数
inline void fluxSolveEuler2D(const std::span<real, 4>& iterVar,
    const std::span<real, 4>& iterFlux,
    const std::array<real, 3>& norm)
{
    //chenyuqing:求解点通量函数
    real r = iterVar[0], u = iterVar[1], v = iterVar[2], p = iterVar[3];
    real Vn = (u * norm[0] + v * norm[1]);
    // real H=GAMMA/(GAMMA-1)*p/r+(u*u+v*v)/2;
    iterFlux[0] = r * Vn;
    iterFlux[1] = r * u * Vn + norm[0] * p;
    iterFlux[2] = r * v * Vn + norm[1] * p;
    //chenyuqing：修改求解点通量函数表达式3
    //iterFlux[2] = r * Vn * v + norm[1] * p;
    iterFlux[3] = ((u * u + v * v) / 2 * r + GAMMA / (GAMMA - 1) * p) * Vn;
    //chenyuqing：修改求解点通量函数表达式4
    //iterFlux[3] = ( GAMMA / (GAMMA - 1) * p + u * u / 2 * r + v * v / 2 * r ) * Vn;
}