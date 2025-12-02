#include "SpaceDis.hpp"
#include "interScheme.hpp"

std::vector<real> SpaceDis::reconR(int i)
{
    std::vector<real> res(nPrim);
    if (interMethod == MUSCL) {
        for (int ivar = 0; ivar < nPrim; ivar++) {
            res[ivar] = musclInterpolation(at(i + 1, ivar), at(i, ivar), at(i - 1, ivar));
        }
        return res;
    } else {
        for (int ivar = 0; ivar < nPrim; ivar++) {
            res[ivar] = at(i, ivar);
        }
        return res;
    }
}

std::vector<real> SpaceDis::reconL(int i)
{
    std::vector<real> res(nPrim);
    if (interMethod == MUSCL) {
        for (int ivar = 0; ivar < nPrim; ivar++) {
            res[ivar] = musclInterpolation(at(i - 2, ivar), at(i - 1, ivar), at(i, ivar));
        }
        return res;
    } else {
        for (int ivar = 0; ivar < nPrim; ivar++) {
            res[ivar] = at(i - 1, ivar);
        }
        return res;
    }
}

std::vector<real> SpaceDis::reconLChar1D(int i)
{
    assert(info->dim == 1);
    assert(info->eqType == EULER);
    enum { R,
        U,
        P };
    real rRef = at(i - 1, R);
    real pRef = at(i - 1, P);
    real cRef = std::sqrt(GAMMA * (pRef / rRef));
    std::array<real, 5> q1, q2, q3;
    for (int j = i - 3; j < i + 2; j++) {
        int iLocal = j - i + 3;
        q1[iLocal] = at(j, R) - at(j, P) / (cRef * cRef);
        // q2[iLocal]=(at(j,P)+at(j,U)*cRef*rRef)/2;
        // q3[iLocal]=(at(j,P)-at(j,U)*cRef*rRef)/2;
        q2[iLocal] = at(j, P);
        q3[iLocal] = at(j, U);
    }
    auto start = std::chrono::high_resolution_clock::now();
    real Q1 = (*this->inter5Positive)(q1);
    real Q2 = (*this->inter5)(q2);
    real Q3 = (*this->inter5)(q3);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;
    // return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
    return { Q1 + (Q2) / (cRef * cRef), Q3, Q2 };
}
std::vector<real> SpaceDis::reconRChar1D(int i)
{
    enum { R,
        U,
        P };
    real rRef = at(i - 1, R);
    real pRef = at(i, P);
    real cRef = std::sqrt(GAMMA * (pRef / rRef));
    std::array<real, 5> q1, q2, q3;
    for (int j = i + 2; j > i - 3; j--) {
        int iLocal = i + 2 - j;
        q1[iLocal] = at(j, R) - at(j, P) / (cRef * cRef);
        // q2[iLocal]=(at(j,P)+at(j,U)*cRef*rRef)/2;
        // q3[iLocal]=(at(j,P)-at(j,U)*cRef*rRef)/2;
        q2[iLocal] = at(j, P);
        q3[iLocal] = at(j, U);
    }
    auto start = std::chrono::high_resolution_clock::now();
    real Q1 = (*this->inter5Positive)(q1);
    real Q2 = (*this->inter5)(q2);
    real Q3 = (*this->inter5)(q3);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;
    // return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
    return { Q1 + (Q2) / (cRef * cRef), Q3, Q2 };
}

static std::array<real, 2> minDif(real u1L, real u2L, real u1R, real u2R)
{
    std::array<real, 3> difs = { std::abs(u1L - u1R), std::abs(u1L - u2R), std::abs(u2L - u1R) };
    int index = std::min_element(difs.begin(), difs.end()) - difs.begin();
    if (index == 0)
        return { u1L, u1R };
    if (index == 1)
        return { u1L, u2R };
    if (index == 2)
        return { u2L, u1R };
    std::cout << "spRecon Error: minDif index out\n";
    return { 0, 0 };
}

static std::array<real, 2> minDif(std::array<real, 3> uL, std::array<real, 3> uR)
{
    std::array<real, 2> difs = { std::abs(uL[0] - uR[0]), std::abs(uL[1] - uR[1]) };
    int index = std::min_element(difs.begin(), difs.end()) - difs.begin();
    if (index == 0)
        return { uL[0], uR[0] };
    else if (index == 1)
        return { uL[1], uR[1] };
    std::cout << "spRecon Error: minDif index out\n";
    return { 0, 0 };
}

static std::array<real, 2> minDif(std::array<real, 2> uL, std::array<real, 2> uR)
{
    // std::array<real,4> difs={std::abs(uL[0]-uR[0]),std::abs(uL[1]-uR[1]),std::abs(uL[0]-uR[1]),std::abs(uL[1]-uR[0])};
    std::array<real, 2> difs = { std::abs(uL[0] - uR[0]), std::abs(uL[1] - uR[1]) };
    int index = std::min_element(difs.begin(), difs.end()) - difs.begin();
    if (index == 0)
        return { uL[0], uR[0] };
    else if (index == 1)
        return { uL[1], uR[1] };
    else if (index == 2)
        return { uL[0], uR[1] };
    else if (index == 3)
        return { uL[1], uR[0] };
    std::cout << "spRecon Error: minDif index out\n";
    return { 0, 0 };
}
// static std::array<real,2> minDif(std::array<real,3> uL,std::array<real,3> uR)
// {
//     std::array<real,4> difs={std::abs(uL[0]-uR[0]),std::abs(uL[1]-uR[1]),std::abs(uL[0]-uR[1]),std::abs(uL[1]-uR[0])};
//     int index=std::min_element(difs.begin(),difs.end())-difs.begin();
//          if(index==0) return{uL[0],uR[0]};
//     else if(index==1)
//     return{uL[1],uR[1]};
//     else if(index==2)
//     return{uL[0],uR[1]};
//     else if(index==3)
//     return{uL[1],uR[0]};
//     std::cout<<"spRecon Error: minDif index out\n";
//     return{0,0};
// }
static std::array<real, 2> minDif2(std::array<real, 3> uL, std::array<real, 3> uR)
{
    std::array<real, 3> difs = { std::abs(uL[0] - uR[0]), std::abs(uL[0] - uR[1]), std::abs(uL[1] - uR[0]) };
    int index = std::min_element(difs.begin(), difs.end()) - difs.begin();
    if (index == 0)
        return { uL[0], uR[0] };
    else if (index == 1)
        return { uL[0], uR[1] };
    else if (index == 2)
        return { uL[1], uR[0] };
    std::cout << "spRecon Error: minDif index out\n";
    return { 0, 0 };
}
std::vector<real> SpaceDis::recon1DFaceCenter(int i)
{
    std::array<real, 3> primL, primR;
    memcpy(&primL[0], &at(i - 1, 0), nVar * sizeof(real));
    memcpy(&primR[0], &at(i, 0), nVar * sizeof(real));
    eigensystemEuler1D eig = eigensystemEuler1D(primL, primR);
    std::array<real, 5> q1L, q2L, q3L, q1R, q2R, q3R;
    for (int j = i - 3; j < i + 3; j++) {
        enum { R,
            U,
            P };
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, P) });

        int iLocal = j - i + 3;
        if (iLocal < 5) {
            q1L[iLocal] = charTemp[0];
            q2L[iLocal] = charTemp[1];
            q3L[iLocal] = charTemp[2];
        }

        iLocal = i + 2 - j;
        if (iLocal < 5) {
            q1R[iLocal] = charTemp[0];
            q2R[iLocal] = charTemp[1];
            q3R[iLocal] = charTemp[2];
        }
    }
    // auto start = std::chrono::steady_clock::now();
    auto Q1LL = inter5(q1L);
    auto Q1RR = inter5(q1R);
    auto Q2LL = inter5(q2L);
    auto Q2RR = inter5(q2R);
    auto Q3LL = inter5(q3L);
    auto Q3RR = inter5(q3R);

    // auto Q1LLTHINC=THINC1(q1L[1],q1L[2],q1L[3]);
    // auto Q1RRTHINC=THINC1(q1R[1],q1R[2],q1R[3]);
    // auto Q2LLTHINC=THINC1(q2L[1],q2L[2],q2L[3]);
    // auto Q2RRTHINC=THINC1(q2R[1],q2R[2],q2R[3]);
    // auto Q3LLTHINC=THINC1(q3L[1],q3L[2],q3L[3]);
    // auto Q3RRTHINC=THINC1(q3R[1],q3R[2],q3R[3]);

    // std::array<real,2> Q1,Q2,Q3;
    // Q1=minDif((std::array<real,2>){Q1LL,Q1LLTHINC},(std::array<real,2>){Q1RR,Q1RRTHINC});
    // Q2=minDif((std::array<real,2>){Q2LL,Q2LLTHINC},(std::array<real,2>){Q2RR,Q2RRTHINC});
    // Q3=minDif((std::array<real,2>){Q3LL,Q3LLTHINC},(std::array<real,2>){Q3RR,Q3RRTHINC});
    // Q1LL=Q1[0];Q1RR=Q1[1];
    // Q2LL=Q2[0];Q2RR=Q2[1];
    // Q3LL=Q3[0];Q3RR=Q3[1];
    // auto stop = std::chrono::steady_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    // timep+=duration;

    auto resTempL = eig.charToPrim({ Q1LL, Q2LL, Q3LL });
    auto resTempR = eig.charToPrim({ Q1RR, Q2RR, Q3RR });
    return { resTempL[0], resTempL[1], resTempL[2], resTempR[0], resTempR[1], resTempR[2] };
}

std::vector<real> SpaceDis::recon2DFaceCenter(int i)
{
    std::array<real, 4> primL, primR;
    memcpy(&primL[0], &at(i - 1, 0), nVar * sizeof(real));
    memcpy(&primR[0], &at(i, 0), nVar * sizeof(real));
    eigensystemEuler2D eig = eigensystemEuler2D(primL, primR, norm);
    std::array<real, 5> q1L, q2L, q3L, q4L, q1R, q2R, q3R, q4R;
    for (int j = i - 3; j < i + 3; j++) {
        enum { R,
            U,
            V,
            P };
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, V), at(j, P) });

        int iLocal = j - i + 3;
        if (iLocal < 5) {
            q1L[iLocal] = charTemp[0];
            q2L[iLocal] = charTemp[1];
            q3L[iLocal] = charTemp[2];
            q4L[iLocal] = charTemp[3];
        }

        iLocal = i + 2 - j;
        if (iLocal < 5) {
            q1R[iLocal] = charTemp[0];
            q2R[iLocal] = charTemp[1];
            q3R[iLocal] = charTemp[2];
            q4R[iLocal] = charTemp[3];
        }
    }
    bool flagL1, flagR1, flagL2, flagR2, flagL3, flagR3, flagL4, flagR4;
    auto start = std::chrono::steady_clock::now();
    auto Q1LL = inter5(q1L);
    auto Q1RR = inter5(q1R);
    auto Q2LL = inter5(q2L);
    auto Q2RR = inter5(q2R);
    auto Q3LL = inter5(q3L);
    auto Q3RR = inter5(q3R);
    auto Q4LL = inter5(q4L);
    auto Q4RR = inter5(q4R);
    auto stop = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;

    auto resTempL = eig.charToPrim({ Q1LL, Q2LL, Q3LL, Q4LL });
    auto resTempR = eig.charToPrim({ Q1RR, Q2RR, Q3RR, Q4RR });
    return { resTempL[0], resTempL[1], resTempL[2], resTempL[3], resTempR[0], resTempR[1], resTempR[2], resTempR[3] };
}
std::vector<real> SpaceDis::reconLChar2D(int i)
{
    assert(info->dim == 2);
    assert(info->eqType == EULER);
    enum { R,
        U,
        V,
        P };
    eigensystemEuler2D eig = eigensystemEuler2D({ at(i - 1, R), at(i - 1, U), at(i - 1, V), at(i - 1, P) }, norm);
    std::array<real, 5> q1, q2, q3, q4;
    for (int j = i - 3; j < i + 2; j++) {
        int iLocal = j - i + 3;
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, V), at(j, P) });
        q1[iLocal] = charTemp[0];
        q2[iLocal] = charTemp[1];
        q3[iLocal] = charTemp[2];
        q4[iLocal] = charTemp[3];
    }
    auto start = std::chrono::high_resolution_clock::now();
    real Q1 = (*this->inter5)(q1);
    real Q2 = (*this->inter5Positive)(q2);
    real Q3 = (*this->inter5)(q3);
    real Q4 = (*this->inter5)(q4);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;
    auto resTemp = eig.charToPrim({ Q1, Q2, Q3, Q4 });
    return { resTemp[0], resTemp[1], resTemp[2], resTemp[3] };
}

std::vector<real> SpaceDis::reconRChar2D(int i)
{
    assert(info->dim == 2);

    enum { R,
        U,
        V,
        P };
    eigensystemEuler2D eig = eigensystemEuler2D({ at(i, R), at(i, U), at(i, V), at(i, P) }, norm);
    std::array<real, 5> q1, q2, q3, q4;
    for (int j = i + 2; j > i - 3; j--) {
        int iLocal = i + 2 - j;
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, V), at(j, P) });
        q1[iLocal] = charTemp[0];
        q2[iLocal] = charTemp[1];
        q3[iLocal] = charTemp[2];
        q4[iLocal] = charTemp[3];
    }
    auto start = std::chrono::high_resolution_clock::now();
    real Q1 = (*this->inter5)(q1);
    real Q2 = (*this->inter5Positive)(q2);
    real Q3 = (*this->inter5)(q3);
    real Q4 = (*this->inter5)(q4);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;
    auto resTemp = eig.charToPrim({ Q1, Q2, Q3, Q4 });
    return { resTemp[0], resTemp[1], resTemp[2], resTemp[3] };
}

// std::vector<real> SpaceDis::reconLChar2D(int i)
// {
//     assert(info->dim==2);
//     assert(info->eqType==EULER);
//     enum{R,U,V,P};
//     double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
//     std::array<real,2> H;
//     real rl=at(i-1,R),ul=at(i-1,U),vl=at(i-1,V),pl=at(i-1,P);
//     real rr=at(i,R),ur=at(i,U),vr=at(i,V),pr=at(i,P);
//     enum{LL,RR};
//     H[LL]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
//     H[RR]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);
//     coef1=std::sqrt(rl);
//     coef2=std::sqrt(rr);
//     real divisor=1.0/(std::sqrt(rl)+std::sqrt(rr));
//     real rRef=std::sqrt(rl*rr);
//     uAvg=(coef1*ul+coef2*ur)*divisor;
//     vAvg=(coef1*vl+coef2*vr)*divisor;
//     HAvg=(coef1*H[LL]+coef2*H[RR])*divisor;
//     q_2Avg=(uAvg*uAvg+vAvg*vAvg)*0.5;
//     real cRef=std::sqrt((GAMMA-1)*(HAvg-q_2Avg));
//     std::array<real,5> q1,q2,q3,q4;
//     for(int j=i-3;j<i+2;j++)
//     {
//         int iLocal=j-i+3;
//         real Vn=(norm[1]>norm[0])?at(j,V):at(j,U);
//         q1[iLocal]=(norm[1]>norm[0])?at(j,U):at(j,V);
//         q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
//         q3[iLocal]=Vn+at(j,P)/(cRef*rRef);
//         q4[iLocal]=-Vn+at(j,P)/(cRef*rRef);
//     }
//     auto start = std::chrono::high_resolution_clock::now();
//     real Q1=(*this->inter5)(q1);
//     real Q2=(*this->inter5Positive)(q2);
//     real Q3=(*this->inter5)(q3);
//     real Q4=(*this->inter5)(q4);
//     auto stop = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
//     timep+=duration;
//     std::vector<real> res={Q2+(Q3+Q4)*0.5*rRef/cRef
//             ,(norm[1]>norm[0]? Q1 : (Q3-Q4)*0.5 )
//             ,(norm[1]>norm[0]? (Q3-Q4)*0.5  : Q1)
//             ,(Q3+Q4)*(cRef*rRef)*0.5};
//     return res;
// }

// std::vector<real> SpaceDis::reconRChar2D(int i)
// {
//     assert(info->dim==2);
//     enum{R,U,V,P};

//     double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
//     std::array<real,2> H;
//     real rl=at(i-1,R),ul=at(i-1,U),vl=at(i-1,V),pl=at(i-1,P);
//     real rr=at(i,R),ur=at(i,U),vr=at(i,V),pr=at(i,P);
//     enum{LL,RR};
//     H[LL]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
//     H[RR]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);
//     coef1=std::sqrt(rl);
//     coef2=std::sqrt(rr);
//     real divisor=1.0/(std::sqrt(rl)+std::sqrt(rr));
//     real rRef=std::sqrt(rl*rr);
//     uAvg=(coef1*ul+coef2*ur)*divisor;
//     vAvg=(coef1*vl+coef2*vr)*divisor;
//     HAvg=(coef1*H[LL]+coef2*H[RR])*divisor;
//     q_2Avg=(uAvg*uAvg+vAvg*vAvg)*0.5;
//     real cRef=std::sqrt((GAMMA-1)*(HAvg-q_2Avg));

//     std::array<real,5> q1,q2,q3,q4;
//     for(int j=i+2;j>i-3;j--)
//     {
//         int iLocal=i+2-j;
//         real Vn=(norm[1]>norm[0])?at(j,V):at(j,U);
//         q1[iLocal]=(norm[1]>norm[0])?at(j,U):at(j,V);
//         q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
//         q3[iLocal]=Vn+at(j,P)/(cRef*rRef);
//         q4[iLocal]=-Vn+at(j,P)/(cRef*rRef);
//     }
//     auto start = std::chrono::high_resolution_clock::now();
//     real Q1=(*this->inter5)(q1);
//     real Q2=(*this->inter5Positive)(q2);
//     real Q3=(*this->inter5)(q3);
//     real Q4=(*this->inter5)(q4);
//     auto stop = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
//     timep+=duration;
//     std::vector<real> res={Q2+(Q3+Q4)*0.5*rRef/cRef
//             ,(norm[1]>norm[0]? Q1 : (Q3-Q4)*0.5 )
//             ,(norm[1]>norm[0]? (Q3-Q4)*0.5  : Q1)
//             ,(Q3+Q4)*(cRef*rRef)*0.5};
//     return res;
// }

std::vector<real> SpaceDis::reconRprim(int i)
{
    std::vector<real> res(nPrim);
    for (int ivar = 0; ivar < nPrim; ivar++) {
        std::array iprims = { at(i + 2, ivar), at(i + 1, ivar), at(i, ivar), at(i - 1, ivar), at(i - 2, ivar) };
        res[ivar] = (*this->inter5)(iprims);
    }
    return res;
}

std::vector<real> SpaceDis::reconLprim(int i)
{
    std::vector<real> res(nPrim);
    for (int ivar = 0; ivar < nPrim; ivar++) {
        std::array iprims = { at(i - 3, ivar), at(i - 2, ivar), at(i - 1, ivar), at(i, ivar), at(i + 1, ivar) };
        res[ivar] = (*this->inter5)(iprims);
    }
    return res;
}
