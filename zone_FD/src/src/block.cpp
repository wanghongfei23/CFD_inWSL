#include"block.hpp"

/**
 * @brief 输出网格为CGNS格式文件的实现
 * 
 * 将当前网格信息写入CGNS格式文件，便于后续可视化或与其他软件交互。
 * 文件名为"grid_c.cgns"，使用Structured网格类型。
 */
void Block::outputCgns()
{
   cgsize_t isize[3][dim];
   int ni,nj,nk,i,j,k;
   int index_file,icelldim,iphysdim,index_base;
   int index_zone,index_coord;
   char basename[33],zonename[33];

/* create gridpoints for simple example: */
   std::vector<real> x;
   x.reserve(nVer);

/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
/* open CGNS file for write */
   if (cg_open("grid_c.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
/* create base (user can give any name) */
   strcpy(basename,"Base");
   icelldim=dim;
   iphysdim=dim;
   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
/* define zone name (user can give any name) */
   strcpy(zonename,"Zone  1");
   for (int idim = 0; idim < dim; idim++)
   {
        /* vertex size */
        isize[0][idim]=iMax[idim];
        /* cell size */
        isize[1][idim]=isize[0][idim]-1;
        /* boundary vertex size (always zero for structured grids) */
        isize[2][idim]=0;
   }
   

/* create zone */
   cg_zone_write(index_file,index_base,zonename,*isize,CGNS_ENUMV(Structured),&index_zone);
/* write grid coordinates (user must use SIDS-standard names here) */
    for (int idim = 0; idim < dim; idim++)
    {
        std::string name=(idim==0)? "CoordinateX":(idim==1)? "CoordinateY":"CoordinateZ";
        for (int i=0 ; i<nVer ; i++ ) x[i]=coorVer(i,idim);
        cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),
        name.c_str(),
        x.data(),&index_coord);
    }
    
   
/* close CGNS file */
   cg_close(index_file);
}

/**
 * @brief 获取单元中心坐标的实现
 * 
 * 通过索引获取指定单元中心的坐标值
 * @param ic 单元索引
 * @param idim 坐标维度索引（0:x, 1:y, 2:z）
 * @return 指定单元在指定维度上的坐标值
 */
real Block::operator()(int ic,int idim)
{
    return coorCel(ic,idim);
}

/**
 * @brief 获取内部单元最大索引数组的实现
 * @return 包含各维度内部单元最大索引的数组
 */
std::array<int,3> Block::getICMax()
{
   return icMax;
}

/**
 * @brief 获取网格最大索引数组的实现
 * @return 包含各维度网格最大索引的数组
 */
std::array<int,3> Block::getIMax()
{
   return iMax;
}

/**
 * @brief 获取网格维度的实现
 * @return 网格维度（1、2或3维）
 */
int Block::getDim()
{
   return dim;
}

/**
 * @brief 获取单元中心坐标数据的实现
 * @param idim 坐标维度索引（0:x, 1:y, 2:z）
 * @return 包含所有单元在指定维度上坐标值的向量
 */
std::vector<real> Block::getCellCoor(int idim)
{
   assert(idim<=dim);
   std::vector<real> res;
   res.reserve(nCel);
   for (int i = 0; i < nCel; i++)
   {
      res.push_back(coorCel(i,idim));
   }
   return res;
   
}

/**
 * @brief 获取顶点坐标数据的实现
 * @param idim 坐标维度索引（0:x, 1:y, 2:z）
 * @return 包含所有顶点在指定维度上坐标值的向量
 */
std::vector<real> Block::getVertexCoor(int idim)
{
   std::vector<real> res;
   res.reserve(nVer);
   int nVerDim=nVer;
   nVerDim/=( dim == 3 ? 1:( dim==2 ? 2 : 4));
   for (int i = 0; i < nVerDim; i++)
   {
      res.push_back(coorVer(i,idim));
   }
   return res;
   
}
std::vector<real> Block::getCellInterval(int i)
{
   std::vector<real> res;
   res.reserve(nVer);
   for (int idim = 0; idim < dim; idim++)
   {
      res.push_back(intervalCel(i,idim));
   }
   return res;
}

real Block::getMinDh(int i)
{
   real res=coorVer(i+1,0)-coorVer(i,0);
   if(dim>=2) res=std::min(coorVer(i+iMax[0],1)-coorVer(i,1),res);
   if(dim>=3) res=std::min(coorVer(i+iMax[0]*iMax[1],2)-coorVer(i,2),res);
   return res;
}