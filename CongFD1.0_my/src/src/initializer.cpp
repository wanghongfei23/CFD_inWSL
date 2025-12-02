#include"initializer.hpp"

void Initializer::solInit(Block* grid,Data* sol)
{
    std::vector<real> tempsol;
    switch (info->eqType)
    {
        /*case begin*/
    case LINEARCONV1D:
        if (grid->dim!=1)
        { std::cout<<"initialize: dim error \n";return;}
        switch (info->nCase)
        {
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                tempsol.push_back(1+sin(x));
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 1:
            //ADR
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                tempsol.push_back(1+sin(x));
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        
        default:
            break;
        }
        break;
    /*case end*/
    
    /*case begin*/
    case BURGERS1D:
        if (grid->dim!=1)
        { std::cout<<"initialize: dim error \n";return;}
        switch (info->nCase)
        {
        /*case 0 begin*/
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                tempsol.push_back(-sin(M_PI*x));
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        /*case 0 end*/
        
        default:
            break;
        }
        break;
    /*case end*/

    case EULER:
        /*case begin*/
        if (grid->dim==1)
        switch (info->nCase)
        {
        /*case 0 begin*/
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=1.4;
                if (x<0)
                {
                    tempsol.push_back(1);
                    tempsol.push_back(0);
                    tempsol.push_back(1.0/(gamma-1)*1);
                }
                else
                {
                    tempsol.push_back(0.125);
                    tempsol.push_back(0);
                    tempsol.push_back(1.0/(gamma-1)*0.1);
                }
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

            case 1:
            //shu-osher
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                if (x<=0.5)
                {
                    tempsol.push_back(3.857143);
                    tempsol.push_back(3.857143*2.629369);
                    tempsol.push_back(1.0/(gamma-1.0)*(10.0+1.0/3.0)
                                     +3.857143*2.629369*2.629369/2);
                }
                else
                {
                    tempsol.push_back(1.0+0.2*sin(5.0*x));
                    tempsol.push_back(0.0);
                    tempsol.push_back(1.0/(gamma-1.0)*1.0);
                }
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

            case 2:
            //Lax
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r,u,p;
                if (x<=0)
                {
                    r=0.445;u=0.698;p=3.528;
                }
                else
                {
                    r=0.5;u=0;p=0.571;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1)*p+r*u*u/2);
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

            case 3:
            //sedov
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                
                real gamma=GAMMA;
                real r,u,E;
                if (std::abs(x)<1e-10)
                {
                    real dx=(*grid)(i,0)-(*grid)(i-1,0);
                    r=1;u=0;E=3200000.0/dx;
                }
                else
                {
                    r=1;u=0;E=1e-12;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*E);
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

            case 4:
            //Woodward-Colella
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r,u,p;
                if (x<0.1)
                {
                    r=1;u=0;p=1000;
                }
                else if (x>=0.9)
                {
                    r=1;u=0;p=100;
                }
                else
                {
                    r=1;u=0;p=0.01;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1)*p+r*u*u/2);
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

        case 5:
            //双稀疏波
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r,u,p;
                if (x<0)
                {
                    r=1;u=-2;p=0.4;
                }
                else
                {
                    r=1;u=2;p=0.4;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1)*p+r*u*u/2);
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        /*case 0 end*/
        default:
            break;
        }
        else if (grid->dim==2)
        switch (info->nCase)
        {
        case 0:
            //2D Riemann Problem;
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                if (x>0.3)
                {
                    if (y>0.3)
                    {
                        tempsol.push_back(1.5);
                        tempsol.push_back(0);
                        tempsol.push_back(0);
                        tempsol.push_back(1.0/(gamma-1)*1.5);
                    }
                    else
                    {
                        tempsol.push_back(0.5323);
                        tempsol.push_back(0);
                        tempsol.push_back(0.5323*1.206);
                        tempsol.push_back(1.0/(gamma-1)*0.3+1.206*1.206/2*0.5323);
                    }
                }
                else
                {
                    if (y>0.3)
                    {
                        tempsol.push_back(0.5323);
                        tempsol.push_back(0.5323*1.206);
                        tempsol.push_back(0);
                        tempsol.push_back(1.0/(gamma-1)*0.3+1.206*1.206/2*0.5323);
                    }
                    else
                    {
                        tempsol.push_back(0.138);
                        tempsol.push_back(0.138*1.206);
                        tempsol.push_back(0.138*1.206);
                        tempsol.push_back(1.0/(gamma-1)*0.029+1.206*1.206*0.138);
                    }
                }
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;
        case 1:
            //2D Riemann Problem vortex;
            //x [-0.5,0.5]
            //y [-0.5,0.5]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (x>0.0)
                {
                    if (y>0.0)
                    {r=1.0;u=0.75;v=-0.5;p=1.0;}
                    else
                    {r=3.0;u=-0.75;v=-0.5;p=1.0;}
                }
                else
                {
                    if (y>0)
                    {r=2.0;u=0.75;v=0.5;p=1.0;}
                    else
                    {r=1.0;u=-0.75;v=0.5;p=1.0;}
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

        case 2:
            //Implosion problem;
            //x [-0.3,0.3]
            //y [-0.3,0.3]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real eps=1e-10;
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (std::abs(x)+std::abs(y)<0.15-eps)
                {
                    r=0.125;
                    u=0.0;
                    v=0.0;
                    p=0.14;
                }
                else
                {
                    r=1.0;
                    u=0.0;
                    v=0.0;
                    p=1.0;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

            
        case 3:
            //R-T instability;
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if(y<=0.5){
                    r=2.0;u=0;p=2.0*y+1.0;
                    real c=std::sqrt(GAMMA*p/r);
                    v=-0.025*c*cos(8.0*M_PI*(x<0.125?x:(0.25-x)));
                }
                else{
                    r=1.0;u=0;p=y+3.0/2.0;
                    real c=std::sqrt(GAMMA*p/r);
                    v=-0.025*c*cos(8.0*M_PI*(x<0.125?x:(0.25-x)));
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

        case 4:
            //Double Mach reflection;
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if(y>=std::sqrt(3)*(x-1.0/6.0)){
                    r=8.0;u=8.25*cos(M_PI/6);v=-8.25*sin(M_PI/6);p=116.5;
                }
                else{
                    r=1.4;u=0;v=0;p=1.0;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

        case 5:
            //2D Riemann Problem another 1;
            //x [-0.5,0.5]
            //y [-0.5,0.5]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (x>0.0)
                {
                    if (y>0.0)
                    {r=0.5313;u=0.0;v=0.0;p=0.4;}
                    else
                    {r=1.0;u=0.0;v=0.7276;p=1.0;}
                }
                else
                {
                    if (y>0)
                    {r=1.0;u=0.7276;v=0.0;p=1.0;}
                    else
                    {r=0.8;u=0.0;v=0.0;p=1.0;}
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;
        
        case 6:
        // 目前还是case5
            //2D Riemann Problem 3;
            //x [-0.5,0.5]
            //y [-0.5,0.5]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (x>0.0)
                {
                    if (y>0.0)
                    {r=0.5323;u=0.0;v=0.0;p=0.4;}
                    else
                    {r=0.8;u=0.0;v=0.0;p=1.0;}
                }
                else
                {
                    if (y>0)
                    {r=1.0;u=0.7276;v=0.0;p=1.0;}
                    else
                    {r=1.0;u=0.0;v=0.7276;p=1.0;}
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;
            
        default:
            break;
        }
        
    /*case end*/
    default:
        break;
    }
}


void Initializer::initUniformBlock(Block* block)
{
    if(!block)
    {
        std::cout<<"Initializer error: empty shared_ptr Block\n";
        return;
    }
    auto iMax=info->iMax;
    int dim=info->dim;
    auto icMax=info->icMax();
    int nVer=1,nCel=1;
    for (int i = 0; i < 3; i++)
    {
        nVer*=iMax[i];
        nCel*=icMax[i];
    }

    block->dim=dim;
    block->nVer=nVer;
    block->nCel=nCel;
    block->icMax=icMax;
    block->iMax=iMax;
    block->coorVer.init(nVer,dim);
    block->coorCel.init(nCel,dim);
    block->intervalCel.init(nCel,dim);
    block->inited=true;

    //consth related
    info->interval=0;
    std::array<real,3> inters;
    double cmin=info->calZone[0];
    double cmax=info->calZone[1];
    // double interval=round((cmax-cmin)/(iMax[0]-1)*1e8)/1e8;
    double interval=(cmax-cmin)/(iMax[0]-1);
    info->interval=interval;
    if (std::abs(*std::max_element(inters.begin(),inters.end())-*std::min_element(inters.begin(),inters.end()))>1e-10)
    std::cout<<"Initalize error: interval incorrect\n";
    info->constH=true;

    //for vertex
    for (int idim = 0; idim < dim; idim++)
    {
        double cmin=info->calZone[idim*2];
        double cmax=info->calZone[idim*2+1];
        double interval=info->interval;

        int l,m,n;
        int* onedIndex=((idim == 0 ) ? &l : ( idim == 1 ? &m : &n));
        
        for (l = 0; l < iMax[0]; l++)
        for (m = 0; m < iMax[1]; m++)
        for (n = 0; n < iMax[2]; n++)
        {
            int globalIndex=l+m*iMax[0]+n*iMax[0]*iMax[1];
            block->coorVer(globalIndex,idim)=(cmin+(*onedIndex)*interval);
            //block->coorVer(globalIndex,idim)=cmin+(*onedIndex)*interval;
        }
    }
    
    //for cellcenter
    // for (int idim = 0; idim < dim; idim++)
    // {
    //     int l,m,n,iLen=(dim==1?2:dim==2? 4:8);
    //     std::vector<int> index;
    //     index.resize(iLen);
    //     for (l = 0; l < icMax[0]; l++)
    //     for (m = 0; m < icMax[1]; m++)
    //     for (n = 0; n < icMax[2]; n++)
    //     {
    //         int iVerGlobal=l+m*iMax[0]+n*iMax[0]*iMax[1];
    //         int iCelGlobal=l+m*icMax[0]+n*icMax[0]*icMax[1];
    //         double temp=0;
    //         index[0]=iVerGlobal;
    //         index[1]=iVerGlobal+1;
    //         if (dim>=2)
    //         {
    //             index[2]=iVerGlobal+iMax[0];
    //             index[3]=iVerGlobal+iMax[0]+1;
    //         }
    //         if (dim>=3)
    //         {
    //             index[4]=iVerGlobal+iMax[0]*iMax[1];
    //             index[5]=iVerGlobal+1+iMax[0]*iMax[1];
    //             index[6]=iVerGlobal+iMax[0]+iMax[0]*iMax[1];
    //             index[7]=iVerGlobal+iMax[0]+1+iMax[0]*iMax[1];
    //         }
    //         real tempInterval;
    //         for(auto iver:index) temp+=block->coorVer(iver,idim);
    //         block->coorCel(iCelGlobal,idim)=temp/iLen;
    //         //only lower line approx_1imate
    //         if(idim==0) block->intervalCel(iCelGlobal,idim)=block->coorVer(index[1],idim)-block->coorVer(index[0],idim);
    //         if(idim==1) block->intervalCel(iCelGlobal,idim)=block->coorVer(index[2],idim)-block->coorVer(index[0],idim);
    //         if(idim==2) block->intervalCel(iCelGlobal,idim)=block->coorVer(index[4],idim)-block->coorVer(index[0],idim);
    //     }
    // }

    for (int idim = 0; idim < dim; idim++)
    {
        int l,m,n,iLen=(dim==1?2:dim==2? 4:8);
        std::vector<int> index;
        index.resize(iLen);
        for (l = 0; l < icMax[0]; l++)
        for (m = 0; m < icMax[1]; m++)
        for (n = 0; n < icMax[2]; n++)
        {
            int iCelGlobal=l+m*icMax[0]+n*icMax[0]*icMax[1];
            int* onedIndex=((idim == 0 ) ? &l : ( idim == 1 ? &m : &n));
            block->coorCel(iCelGlobal,idim)=(cmin+((*onedIndex)+0.5)*interval);
            //only lower line approx_1imate
            if(idim==0) block->intervalCel(iCelGlobal,idim)=interval;
            if(idim==1) block->intervalCel(iCelGlobal,idim)=interval;
            if(idim==2) block->intervalCel(iCelGlobal,idim)=interval;
        }
    }
}


Initializer::Initializer()
{

}

Initializer::Initializer(Info* info_)
{
    info=info_;
}


void Initializer::initEqution(Equation* eq,Block* block)
{
    if(!eq||!block)
    {
        std::cout<<"Initializer error: empty shared_ptr Equation or Block\n";
        return;
    }
    eq->n=block->icMax[0]*block->icMax[1]*block->icMax[2];
    eq->nPrim=info->nPrim();
    eq->nCons=info->nCons();
    eq->type=info->eqType;
    eq->dim=block->dim;

    
    eq->rhs=new Data(eq->n,eq->nCons);
    eq->cons=new Data(eq->n,eq->nCons);
    eq->prim=new Data(eq->n,eq->nPrim);
    eq->inited=true;
    eq->cons->setvarName(info->getVarNameListCons());
    eq->prim->setvarName(info->getVarNameListPrim());
    eq->rhs->setvarName(info->getVarNameListRhs());

    solInit(block,eq->cons);
}


void Initializer::initBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block)
{
    bnds->iMax=iMax;
    bnds->dim=info->dim;

    int nGhost=info->nGhostCell();
    int nCons=info->nCons();
    int nPrim=info->nPrim();

    std::array<int,3> nBnds  {bnds->iMax[1]*bnds->iMax[2],
                            bnds->iMax[0]*bnds->iMax[2],
                            bnds->iMax[0]*bnds->iMax[1]};
    int nBnd=0;
    for (int i = 0; i < bnds->dim; i++)
    {
        nBnd+=nBnds[i]*2;
    }
    
    

    bnds->oneDBnds.resize(nBnd);
    std::array<int,2> offsets;
    switch (info->eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        {
            for (int i = 0; i < 2; i++)
            {
                bnds->oneDBnds.at(i)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
            }
            offsets=calOffsetInverse(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(0)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            offsets=calOffset(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }
        break;
    case EULER:
        if(eqn->dim==1)
        {
            BndType Xtype=SUPERSONICOUTLET;
            if(info->nCase==4) {Xtype=SYMMETRY1D;}
            bnds->oneDBnds.at(0)=std::make_shared<OneDBnd>(nGhost,nPrim,Xtype);
            offsets=calOffset(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(0)->setUpdate(eqn->prim,offsets[0],offsets[1]);


            bnds->oneDBnds.at(1)=std::make_shared<OneDBnd>(nGhost,nPrim,Xtype);
            offsets=calOffsetInverse(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }
        else if (eqn->dim==2)
        {
            BndType Xtype=SUPERSONICOUTLET,Ytype=SUPERSONICOUTLET;
            if(info->nCase==2) {Xtype=SYMMETRYX;Ytype=SYMMETRYY;}
            if(info->nCase==3) {Xtype=SYMMETRYX;Ytype=DIRICLET;}
            if(info->nCase==4) 
            {
                initDoubleMachBnds(bnds,eqn,iMax,block);
                break;
            }

            for (int i = 0; i < iMax[1]; i++)
            {
                bnds->oneDBnds.at(2*i)=std::make_shared<OneDBnd>(nGhost,nPrim,Xtype);
                offsets=calOffset(1,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                bnds->oneDBnds.at(2*i+1)=std::make_shared<OneDBnd>(nGhost,nPrim,Xtype);
                offsets=calOffsetInverse(1,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
            }

            for (int i = 0; i < iMax[0]; i++)
            {
                bnds->oneDBnds.at(2*i+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,Ytype);
                offsets=calOffset(2,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                if(info->nCase==3) //For R-T instability
                {
                    std::array<real,4> dirVar={2.0,0,0,1.0};
                    std::vector<real> dirVars(nGhost*nPrim);
                    for(int j=0;j<nGhost;j++) {
                        dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                        dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
                    }
                    bnds->oneDBnds.at(2*i+iMax[1]*2)->setValue(dirVars);
                }

                bnds->oneDBnds.at(2*i+1+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,Ytype);
                offsets=calOffsetInverse(2,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                if(info->nCase==3) //For R-T instability
                {
                    std::array<real,4> dirVar={1.0,0,0,2.5};
                    std::vector<real> dirVars(nGhost*nPrim);
                    for(int j=0;j<nGhost;j++) {
                        dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                        dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
                    }
                    bnds->oneDBnds.at(2*i+iMax[1]*2+1)->setValue(dirVars);
                }
            }
            
        }
        

        break;
    
    default:
        break;
    }
    
}

void Initializer::initDoubleMachBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block)
{
        std::array<int,2> offsets;
        int nGhost=info->nGhostCell();
        int nCons=info->nCons();
        int nPrim=info->nPrim();
        for (int i = 0; i < iMax[1]; i++)
        {
            offsets=calOffset(1,i,0,bnds->iMax);
            real x=block->coorCel(offsets[0],0);
            bnds->oneDBnds.at(2*i)=std::make_shared<OneDBnd>(nGhost,nPrim,DIRICLET);
            bnds->oneDBnds.at(2*i)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            //set values for dirichlet boundary
            std::array<real,4> dirVar={8.0,8.25*cos(M_PI/6),-8.25*sin(M_PI/6),116.5};
            std::vector<real> dirVars(nGhost*nPrim);
            for(int j=0;j<nGhost;j++) {
                dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
            }
            bnds->oneDBnds.at(2*i)->setValue(dirVars);

            bnds->oneDBnds.at(2*i+1)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
            offsets=calOffsetInverse(1,i,0,bnds->iMax);
            bnds->oneDBnds.at(2*i+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }

        for (int i = 0; i < iMax[0]; i++)
        {
            offsets=calOffset(2,i,0,bnds->iMax);
            real x=block->coorCel(offsets[0],0);
            BndType Ytype=SYMMETRYY;
            if(x<1.0/6.0) Ytype=DIRICLET;
            bnds->oneDBnds.at(2*i+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,Ytype);
            bnds->oneDBnds.at(2*i+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            if(Ytype==DIRICLET)
            {
                std::array<real,4> dirVar={8.0,8.25*cos(M_PI/6),-8.25*sin(M_PI/6),116.5};
                std::vector<real> dirVars(nGhost*nPrim);
                for(int j=0;j<nGhost;j++) {
                    dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                    dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
                }
                bnds->oneDBnds.at(2*i+iMax[1]*2)->setValue(dirVars);
            }

            bnds->oneDBnds.at(2*i+1+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,DoubleMachUp);
            offsets=calOffsetInverse(2,i,0,bnds->iMax);
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            std::array<real,3> coor={block->coorCel(offsets[0],0),
                                    block->coorCel(offsets[0],1),
                                    block->coorCel(offsets[0],2)};
            std::array<real,3> dh={block->coorCel(offsets[0],0)-block->coorCel(offsets[0]+offsets[1],0),
                                   block->coorCel(offsets[0],1)-block->coorCel(offsets[0]+offsets[1],1),
                                   block->coorCel(offsets[0],2)-block->coorCel(offsets[0]+offsets[1],2)};
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setInfo(info);
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setCoor(coor,dh);
        }
}

void Initializer::initSpDistributor(SpDistributor* spDis,Equation* eqn
                                    ,Block* block,Bnds* bnds)
{
    
    spDis->nCons=info->nCons();
    spDis->nPrim=info->nPrim();
    spDis->iMax=block->icMax;
    spDis->prim=eqn->prim;
    spDis->cons=eqn->cons;
    spDis->dim=info->dim;
    spDis->bnds=bnds;
    spDis->rhs=eqn->rhs;
    spDis->info=info;
    

}