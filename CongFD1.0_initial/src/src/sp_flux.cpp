#include "SpaceDis.hpp"
#include "fluxScheme.hpp"
#include "interScheme.hpp"

void SpaceDis::calFluxConv(int i)
{
    /*for u_t + a * u_x == 0*/
    real a = 1.0;
    real ul, ur;
    ul = (this->*reconLMethod)(i).at(0);
    ur = (this->*reconRMethod)(i).at(0);
    /*
    ul=data(i-1,0);
    ur=data(i,0);*/
    /*if (a>0)
    {
        flux(i,0)=a*ul;
    }
    else
    {
        flux(i,0)=a*ur;
    }*/

    fluxAt(i, 0) = 0.5 * (a * ul + a * ur - std::abs(a) * (ur - ul));
}
void SpaceDis::calFluxAccuracyTest(int i)
{
    /*for u_t + a * u_x == 0*/
    real a = 1.0;
    real ul, ur;
    ul = (this->*reconLMethod)(i).at(0);
    fluxAt(i, 0) = ul;
}

void SpaceDis::calFluxBurgers(int i)
{
    /*for u_t + a * u_x == 0*/

    real ul, ur, aLF;
    aLF = data->maxElement(0);

    ul = (this->*reconLMethod)(i).at(0);
    ur = (this->*reconRMethod)(i).at(0);

    // ul=data(i-1,0);
    // ur=data(i,0);
    // L-F flux
    real al, ar;
    al = std::abs(ul);
    ar = std::abs(ur);
    // real a=std::max(ar,ar);
    // real a=(al+ar)/2;
    real a = aLF;
    fluxAt(i, 0) = 0.5 * (ul * ul / 2 + ur * ur / 2 - a * (ur - ul));
}

typedef std::array<real, 2> arr2;

// void SpaceDis::calFluxEuler1D(int i)
// {
//     /*for u_t + a * u_x == 0*/
//     auto WL=(this->*reconLMethod)(i);
//     auto WR=(this->*reconRMethod)(i);
//     enum{R,U,P};
//     if (WL[P]<0||WL[R]<0)
//     {
//         WL[R]=at(i-1,R);
//         WL[U]=at(i-1,U);
//         WL[P]=at(i-1,P);
//         std::cout<<"SpaceDis error: negative pressure/Density \n";
//     }
//     if (WR[P]<0||WR[R]<0)
//     {
//         WR[R]=at(i,R);
//         WR[U]=at(i,U);
//         WR[P]=at(i,P);
//         std::cout<<"SpaceDis error: negative pressure/Density \n";
//     }
//     std::vector<real> iflux=roeFlux1D2(WL[0],WR[0],WL[1],WR[1],WL[2],WR[2]);

//     //auto iflux=HLLCFlux1D(WL[0],WR[0],WL[1],WR[1],WL[2],WR[2]);
//     //std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
//     for (int ivar = 0; ivar < 3; ivar++)
//     {
//         fluxAt(i,ivar)=iflux[ivar];
//     }
// }

void SpaceDis::calFluxEuler1D(int i)
{
    /*for u_t + a * u_x == 0*/
    auto W = this->recon1DFaceCenter(i);
    enum { R,
        U,
        P };
    if (W[0] < 0 || W[2] < 0) {
        W[0] = at(i - 1, R);
        W[1] = at(i - 1, U);
        W[2] = at(i - 1, P);
        std::cout << "SpaceDis error: negative pressure/Density \n";
    }
    if (W[3] < 0 || W[5] < 0) {
        W[3] = at(i, R);
        W[4] = at(i, U);
        W[5] = at(i, P);
        std::cout << "SpaceDis error: negative pressure/Density \n";
    }
    std::vector<real> iflux = roeFlux1D2(W[0], W[3], W[1], W[4], W[2], W[5]);

    // auto iflux=HLLCFlux1D(W[0],W[3],W[1],W[4],W[2],W[5]);
    // std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
    for (int ivar = 0; ivar < 3; ivar++) {
        fluxAt(i, ivar) = iflux[ivar];
    }
}

void SpaceDis::calFluxEuler2D(int i)
{
    /*for u_t + a * u_x == 0*/

    // auto WL=(this->*reconLMethod)(i);
    // auto WR=(this->*reconRMethod)(i);
    auto W = this->recon2DFaceCenter(i);
    // if (WL[3]<0||WL[0]<0)//||isnan(WL[3])||isnan(WL[0]))
    if (W[3] < 0 || W[0] < 0) {
        // int ivar=(WL[3]<0||isnan(WL[3]))? 3:0;
        // std::array<real,5> q={at(i-3,ivar),at(i-2,ivar),at(i-1,ivar),at(i,ivar),at(i+1,ivar)};
        // std::array<real,3> beta;
        // unsigned flag=0;
        // beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
        //     + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);
        // beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
        //     + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

        // beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
        //     + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

        // auto minBeta=std::min_element(beta.begin(),beta.end())-beta.begin();
        // //if (beta[minBeta]<eps) return q[2];
        // real CT=CTi,sumGamma=0;
        // real tau=std::abs(beta[2]-beta[0]),KK=CT*(beta[minBeta]+tau);
        // real eps=1e-8;
        // if (minBeta==0||(beta[0]+tau)*beta[minBeta]>KK*beta[0]+eps)
        // {
        //     flag+=1;
        // }
        // if (minBeta==1||(beta[1]+tau)*beta[minBeta]>KK*beta[1]+eps)
        // {
        //     flag+=2;
        // }
        // if (minBeta==2||(beta[2]+tau)*beta[minBeta]>KK*beta[2]+eps)
        // {
        //     flag+=4;
        // }

        // std::cout<<std::format("SpaceDis error L: negative pressure i={} WLR={} WRR={} WLP={} WRP={} \
        //                         WSR={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} \
        //                         WSU={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} \
        //                         WSV={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} \
        //                         WSP={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} \
        //                         beta={:.4f} {:.4f} {:.4f} flag={}\n"
        //                       ,i,WL[0],WR[0],WL[3],WR[3]
        //                       ,at(i-3,0),at(i-2,0),at(i-1,0),at(i,0),at(i+1,0)
        //                       ,at(i-3,1),at(i-2,1),at(i-1,1),at(i,1),at(i+1,1)
        //                       ,at(i-3,2),at(i-2,2),at(i-1,2),at(i,2),at(i+1,2)
        //                       ,at(i-3,3),at(i-2,3),at(i-1,3),at(i,3),at(i+1,3)
        //                       ,beta[0],beta[1],beta[2],flag);
        std::cout << std::format("SpaceDis error L: negative pressure i={}\n", i);
        W = { at(i - 1, 0), at(i - 1, 1), at(i - 1, 2), at(i - 1, 3), W[4], W[5], W[6], W[7] };
        auto WX = (this->*reconLMethod)(i);
    }
    // if(W[3]<0||W[0]<0)//||isnan(WR[3])||isnan(WR[0]))
    if (W[3 + 4] < 0 || W[0 + 4] < 0) //||isnan(WR[3])||isnan(WR[0]))
    {
        // std::array<real,5> q={at(i+2,3),at(i+1,3),at(i,3),at(i-1,3),at(i-2,3)};
        // std::array<real,3> beta;
        // beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
        //     + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);
        // beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
        //     + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);
        // beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
        //     + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
        // std::cout<<std::format("SpaceDis error R: negative pressure i={} WLR={} WRR={} WLP={} WRP={} WS={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} beta={:.4f} {:.4f} {:.4f}\n"
        //                       ,i,WL[0],WR[0],WL[3],WR[3],at(i-2,3),at(i-1,3),at(i,3),at(i+1,3),at(i+2,3),beta[0],beta[1],beta[2]);
        std::cout << std::format("SpaceDis error R: negative pressure i={}\n", i);
        W = { W[0], W[1], W[2], W[3], at(i, 0), at(i, 1), at(i, 2), at(i, 3) };
        // auto WX=(this->*reconRMethod)(i);
    }

    // auto iflux = roeFlux2DSym(W[0], W[4], W[1], W[5], W[2], W[6], W[3], W[7], norm);
     auto iflux=HLLCFlux2D2(W[0],W[4],W[1],W[5],W[2],W[6],W[3],W[7],norm);
    // std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
    for (int ivar = 0; ivar < 4; ivar++) {
        fluxAt(i, ivar) = iflux[ivar];
    }
}
