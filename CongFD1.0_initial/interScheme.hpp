#define real double
#include <array>
#include <cmath>
#include <algorithm>

real Teno5_Z(std::array<real,5> q)
{
    real eps=1e-40;
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

    real sumbeta=0;
    real C=1,qq=6,tau=std::abs(beta[2]-beta[0]);
    for(int i=0;i<3;i++)
    {
        real tempp=C+tau/(beta[i]+eps);
        tempp*=tempp;
        beta[i]=tempp*tempp*tempp;
        sumbeta+=beta[i];
    }
    real CT=1e-5*sumbeta;
    
    unsigned short flag=0;
    if(beta[0]<CT) flag+=1;
    if(beta[1]<CT) flag+=2;
    if(beta[2]<CT) flag+=4;
    switch (flag)
    {
    case 0:
        /* 1,1,1 */
        return 3.0/128.0*q[0]-5.0/32.0*q[1]+45.0/64.0*q[2]+15.0 /32.0*q[3]-5.0/128.0*q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0/16.0*q[1]+9.0/16.0*q[2]+9.0 /16.0*q[3]-1.0/16.0*q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0/16.0*q[0]-5.0/16.0*q[1]+15.0/16.0*q[2]+5.0/16.0*q[3];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

real Teno5_SZ(std::array<real,5> q)
{
    real eps=1e-40;//1e-10;
    std::array<real,3> beta={
    1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2),

    1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2),

    1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2)};

    unsigned short minBeta=std::min_element(beta.begin(),beta.end())-beta.begin();
    constexpr real CT=0.15704178024750198;        //C_T=1e-5
    //constexpr real CT=0.10699131939336631;      //C_T=1e-6
    //constexpr real CT=0.072892337360747711;     //C_T=1e-7
    //constexpr real CT=0.033833625914958219;     //C_T=1e-9
    constexpr real CT_1=1-CT;
    real tau=std::abs(beta[2]-beta[0]);
    real rr=CT*tau-CT_1*beta[minBeta];
    real ll=tau*beta[minBeta];

    unsigned short flag=0;
    if(ll<rr*beta[0]) flag+=1;
    if(ll<rr*beta[1]) flag+=2;
    if(ll<rr*beta[2]) flag+=4;
    switch (flag)
    {
    case 0:
        /* 1,1,1 */
        return 3.0/128.0*q[0]-5.0/32.0*q[1]+45.0/64.0*q[2]+15.0 /32.0*q[3]-5.0/128.0*q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0/16.0*q[1]+9.0/16.0*q[2]+9.0 /16.0*q[3]-1.0/16.0*q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0/16.0*q[0]-5.0/16.0*q[1]+15.0/16.0*q[2]+5.0/16.0*q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

real weno5_JS(std::array<real,5> q)
{
    real eps=1e-6;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);


    std::array<real,3> u;
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    
    real sumbeta=0,result=0;
    for(int i=0;i<3;i++)
    {
        beta[i]=gamma[i]/pow(eps+beta[i],2);
        sumbeta+=beta[i];
    }
    for(int i=0;i<3;i++) result+=beta[i]*u[i];
    return result/sumbeta;
}

real weno5_Z(std::array<real,5> q)
{
    real eps=1e-40;
    std::array<real,3> gamma={1.0/16.0,5.0/8.0,5.0/16.0};
    std::array<real,3> beta;
    beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
            + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
            + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
            + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

    std::array<real,3> u;
    u[0]= 3.0/8.0*q[0]-5.0/4.0*q[1]+15.0/8.0*q[2];
    u[1]=-1.0/8.0*q[1]+3.0/4.0*q[2]+3.0 /8.0*q[3];
    u[2]= 3.0/8.0*q[2]+3.0/4.0*q[3]-1.0 /8.0*q[4];
    
    real sumbeta=0,result=0;
    real C=1,qq=2,tau=std::abs(beta[2]-beta[0]);
    for(int i=0;i<3;i++)
    {
        beta[i]=gamma[i]*(C+pow(tau/(beta[i]+eps),qq));
        sumbeta+=beta[i];
    }
    for(int i=0;i<3;i++) result+=beta[i]*u[i];
    return result/sumbeta;
}