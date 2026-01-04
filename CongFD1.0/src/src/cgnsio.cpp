#include "cgnsio.hpp"

void CgnsIO::BlockCgnsOutput(Block* block,Info* info)
{
   {
   int dim=block->getDim();
   
   auto iMax=block->getIMax();
   cgsize_t isize[3][dim]={0};
   auto filename=info->filename();
   int index_file,icelldim,iphysdim,index_base;
   int index_zone,index_coord;
   char basename[33],zonename[33];

/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
/* open CGNS file for write */
   if (cg_open(filename.c_str(),CG_MODE_WRITE,&index_file)) cg_error_exit();
/* create base (user can give any name) */
   strcpy(basename,"Base");
   icelldim=dim;
   iphysdim=dim;
   if (cg_base_write(index_file,basename,icelldim,iphysdim,&index_base)) cg_error_exit();
/* define zone name (user can give any name) */
   strcpy(zonename,"Zone  1");
   for (int idim = 0; idim < dim; idim++)
   {
        /* vertex size */
        isize[0][idim]=iMax[idim];
        /* cell size */
        isize[1][idim]=iMax[idim]-1;
        /* boundary vertex size (always zero for structured grids) */
        isize[2][idim]=0;
   }
   

/* create zone */
   if (cg_zone_write(index_file,index_base,zonename,*isize,CGNS_ENUMV(Structured),&index_zone)) cg_error_exit();
/* write grid coordinates (user must use SIDS-standard names here) */
    for (int idim = 0; idim < dim; idim++)
    {
        std::string name=(idim==0)? "CoordinateX":(idim==1)? "CoordinateY":"CoordinateZ";
        auto x=block->getVertexCoor(idim);
        if (cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),
        name.c_str(),
        x.data(),&index_coord)) cg_error_exit();
    }
    
   
/* close CGNS file */
   cg_close(index_file);}
}



void CgnsIO::solCgnsOutput(Data* data,Info* info)
{

    int index_file;
    auto filename=info->filename();
/* WRITE FLOW SOLUTION TO EXISTING CGNS FILE */
/* open CGNS file for modify */
    if (cg_open(filename.c_str(),CG_MODE_MODIFY,&index_file)) cg_error_exit();
/* we know there is only one base (real working code would check!) */
   int index_base,index_zone,index_flow,index_field;
    index_base=1;
/* we know there is only one zone (real working code would check!) */
    index_zone=1;
/* define flow solution node name (user can give any name) */
   std::string solname="flowSolution";
/* create flow solution node */
    if (cg_sol_write(index_file,index_base,index_zone,solname.c_str(),CGNS_ENUMV(CellCenter),&index_flow)) cg_error_exit();
/* write flow solution (user must use SIDS-standard names here) */
   for(int ivar=0;ivar<data->getNVar();ivar++)
   {
      auto varr=data->getIVar(ivar);
      if (cg_field_write(index_file,index_base,index_zone,index_flow,
        CGNS_ENUMV(RealDouble),data->varName[ivar].c_str(),varr.data(),&index_field)) cg_error_exit();
   }

/* close CGNS file */
    if (cg_close(index_file)) cg_error_exit();;
}
