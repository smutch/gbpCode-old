#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>

int write_treenode_list_properties_set_ith(tree_info *trees,int i_write,tree_node_info *current_halo,char *data_name,SID_Datatype *data_type,int *data_i,double *data_d){
   int i_item=0;
   if(i_write==(i_item++)){
      // We always need to check pointers first for the case when we're calling this function
      //    to count entries or to write the header 
      if(data_name!=NULL) sprintf(data_name,"Expansion factor");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_i!=NULL)    *data_d   =trees->a_list[current_halo->snap_tree];
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Redshift");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_i!=NULL)    *data_d   =trees->z_list[current_halo->snap_tree];
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Snapshot");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =trees->snap_list[current_halo->snap_tree];
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"File index");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =current_halo->file_index;
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Halo ID");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_halo_ID(trees,current_halo);
   } 
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Descendant ID");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_halo_ID(trees,current_halo->descendant);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"No. of particles");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_n_particles(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"No. of particles peak");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_n_particles_peak(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"No. of particles inclusive");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_n_particles_inclusive(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"No. of particles inclusive peak");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_n_particles_inclusive_peak(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"No. of particles descendant");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_n_particles(trees,current_halo->descendant);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"No. of particles descendant peak");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_n_particles_peak(trees,current_halo->descendant);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"No. of particles top parent");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_n_particles(trees,current_halo->parent_top);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Substructure rank");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_substructure_rank(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"M_vir [h^{-1} M_sol]");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =fetch_treenode_M_vir(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"M_peak [h^{-1} M_sol]");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =fetch_treenode_M_peak(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Progenitor rank");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =fetch_treenode_progenitor_rank(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"x-position [h^{-1} Mpc]");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =fetch_treenode_x(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"y-position [h^{-1} Mpc]");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =fetch_treenode_y(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"z-position [h^{-1} Mpc]");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =fetch_treenode_z(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Separation from main progenitor [units of main progenitor R_vir]");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =fetch_treenode_delta_r_MP(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Branch age [Gyrs]");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =fetch_treenode_delta_t_leaf(trees,current_halo)/S_PER_GYR;
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Half-mass age [Gyrs]");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =fetch_treenode_delta_t_form(trees,current_halo)/S_PER_GYR;
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Descendant match score");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =(double)fetch_treenode_descendant_score(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Progenitor match score");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =(double)fetch_treenode_progenitor_score(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Descendant f_goodness");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =(double)fetch_treenode_descendant_f_goodness(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Progenitor f_goodness");
      if(data_type!=NULL) *data_type=SID_DOUBLE;
      if(data_d!=NULL)    *data_d   =(double)fetch_treenode_progenitor_f_goodness(trees,current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Fragmented?");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =check_treenode_if_fragmented(current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Central?");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =check_treenode_if_central(current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Satellite?");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =check_treenode_if_satellite(current_halo);
   }
   else if(i_write==(i_item++)){
      if(data_name!=NULL) sprintf(data_name,"Tree case BWS");
      if(data_type!=NULL) *data_type=SID_INT;
      if(data_i!=NULL)    *data_i   =current_halo->tree_case;
   } 
   // 'else' is structured as it is so that you can iteratively call this function (with all pointers==NULL)
   //    until it returns FALSE, to get a count of how many entries it supports.  This allows for easy changes
   //    to the output format of this file.
   else{ 
      if(data_d!=NULL || data_i!=NULL || data_name!=NULL) 
         SID_trap_error("Invalid counter given in write_treenode_list_properties_set_ith()",ERROR_LOGIC);
      else
         return(FALSE);
   }
   return(TRUE);
}

