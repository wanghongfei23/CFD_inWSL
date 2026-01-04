#pragma once
#include<array>
#include<macro.hpp>

enum{X,Y};
class eigensystemEuler2D
{
    public:
    eigensystemEuler2D(){};
    eigensystemEuler2D(const std::array<real,4> & prim,const std::array<real,3> & norm_);
    eigensystemEuler2D(const std::array<real,4> &priml,const std::array<real,4> &primr,const std::array<real,3> & norm_);
    std::array<real,4> primToChar(const std::array<real,4> & prim);
    std::array<real,4> charToPrim(const std::array<real,4> & chars);
    

    private:
    real r,u,v,p,gamma=GAMMA,ek,h,c,Vn;
    bool xOrY=false;
    std::array<real,3> norm;
    std::array<real,4*4> leftEig,rightEig;
};

class eigensystemEuler1D
{
    public:
    eigensystemEuler1D(const std::array<real,3> & prim);
    eigensystemEuler1D(const std::array<real,3> &priml,const std::array<real,3> &primr);
    std::array<real,3> primToChar(const std::array<real,3> & prim);
    std::array<real,3> charToPrim(const std::array<real,3> & chars);
    

    private:
    real r,u,p,gamma=GAMMA,ek,h,c;
    bool xOrY=false;
    std::array<real,3*3> leftEig,rightEig;
};

// std::array<real,4> characteristicDecomposition(std::array<real,4> prim)
// {
//     return{0,0,0,0};
// }
// std::array<real,4*4> leftEigen(const std::array<real,4>& prim,const std::array<real,3>& norm)
// {
//     std::array<real,4*4> res;
    
//     real r=prim[0],u=prim[1],v=prim[2],p=prim[3];
//     real gamma=GAMMA,ek=(u*u+v*v)/2;
//     real h=p/r*gamma/(1-gamma);
//     real Vn=norm[0]*u+norm[1]*v;
//     real c=std::sqrt(gamma*p/r);
//     //first line
//     res[0]=-norm[Y]*u+norm[X]*v;
//     res[1]=norm[Y];
//     res[2]=-norm[X];
//     res[3]=0;

//     //second line
//     res[4]=h-ek;
//     res[5]=u;
//     res[6]=v;
//     res[7]=-1;

//     //third line
//     res[8 ]=(Vn/c+ek/h)/2;
//     res[9 ]=(-norm[X]/c-u/h)/2;
//     res[10]=(-norm[Y]/c-v/h)/2;
//     res[11]=1.0/(2*h);

//     //fourth line
//     res[12]=(-Vn/c+ek/h)/2;
//     res[13]=(norm[X]/c-u/h)/2;
//     res[14]=(norm[Y]/c-v/h)/2;
//     res[15]=1.0/(2*h);

//     return res;
// }

// std::array<real,4*4> rightEigen(const std::array<real,4>& prim,const std::array<real,3>& norm)
// {
//     std::array<real,4*4> res;
//     enum{X,Y};
//     real r=prim[0],u=prim[1],v=prim[2],p=prim[3];
//     real gamma=GAMMA,ek=(u*u+v*v)/2;
//     real h=p/r*gamma/(1-gamma);
//     real Vn=norm[0]*u+norm[1]*v;
//     real c=std::sqrt(gamma*p/r);
//     //first line
//     res[0]=0;
//     res[1]=1/h;
//     res[2]=1;
//     res[3]=1;

//     //second line
//     res[4]=norm[Y];
//     res[5]=u/h;
//     res[6]=u-norm[X]*c;
//     res[7]=u+norm[X]*c;

//     //third line
//     res[8 ]=-norm[X];
//     res[9 ]=v/h;
//     res[10]=v-norm[Y]*c;
//     res[11]=v+norm[Y]*c;

//     //fourth line
//     res[12]=norm[Y]*u-norm[X]*v;
//     res[13]=ek/h;
//     res[14]=h+ek-Vn*c;
//     res[15]=h+ek+Vn*c;
//     return res;
// }