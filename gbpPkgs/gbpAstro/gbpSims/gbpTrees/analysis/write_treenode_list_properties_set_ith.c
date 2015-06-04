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

// N.B.: DON'T FORGET TO ADJUST WRITE_TREENODE_LIST_PROPERTIES_N IF YOU ADD/REMOVE ITEMS FROM THIS LIST!

void write_treenode_list_properties_set_ith(tree_info *trees,tree_node_info *current_halo,int i_write,char *data_name,SID_Datatype *data_type,int *data_i,double *data_d){
   switch(i_write){
      case 0:
         if(data_name!=NULL) sprintf(data_name,"snapshot");
         if(data_type!=NULL) *data_type=SID_INT;
         if(data_i!=NULL)    *data_i   =trees->snap_list[current_halo->snap_tree];
         break;
      case 1:
         if(data_name!=NULL) sprintf(data_name,"catalog index");
         if(data_type!=NULL) *data_type=SID_INT;
         if(data_i!=NULL)    *data_i   =current_halo->file_index;
         break;
      case 2:
         if(data_name!=NULL) sprintf(data_name,"No. of particles");
         if(data_type!=NULL) *data_type=SID_INT;
         if(data_i!=NULL)    *data_i   =fetch_treenode_n_particles(trees,current_halo);
         break;
      case 3:
         if(data_name!=NULL) sprintf(data_name,"M_vir [h^{-1} M_sol]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_Mvir(trees,current_halo);
         break;
      case 4:
         if(data_name!=NULL) sprintf(data_name,"M_peak [h^{-1} M_sol]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_Mpeak(trees,current_halo);
         break;
      case 5:
         if(data_name!=NULL) sprintf(data_name,"M_peak_descendant [h^{-1} M_sol]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_Mpeak(trees,current_halo->descendant);
         break;
      case 6:{
         if(data_name!=NULL) sprintf(data_name,"M_peak_main_progenitor [h^{-1} M_sol]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_Mpeak(trees,current_halo->descendant->progenitor_first);
         break;
      }
      case 7:
         if(data_name!=NULL) sprintf(data_name,"Progenitor rank");
         if(data_type!=NULL) *data_type=SID_INT;
         if(data_i!=NULL)    *data_i   =fetch_treenode_progenitor_rank(trees,current_halo);
         break;
      case 8:
         if(data_name!=NULL) sprintf(data_name,"x-position [h^{-1} Mpc]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_x(trees,current_halo);
         break;
      case 9:
         if(data_name!=NULL) sprintf(data_name,"y-position [h^{-1} Mpc]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_y(trees,current_halo);
         break;
      case 10:
         if(data_name!=NULL) sprintf(data_name,"z-position [h^{-1} Mpc]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_z(trees,current_halo);
         break;
      case 11:
         if(data_name!=NULL) sprintf(data_name,"x-velocity [km/s]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_vx(trees,current_halo);
         break;
      case 12:
         if(data_name!=NULL) sprintf(data_name,"y-velocity [km/s]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_vy(trees,current_halo);
         break;
      case 13:
         if(data_name!=NULL) sprintf(data_name,"z-velocity [km/s]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_vz(trees,current_halo);
         break;
      case 14:
         if(data_name!=NULL) sprintf(data_name,"Separation from main progenitor [units of main progenitor R_vir]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_delta_r_MP(trees,current_halo);
         break;
      case 15:
         if(data_name!=NULL) sprintf(data_name,"Branch age [Gyrs]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_delta_t_leaf(trees,current_halo)/S_PER_GYR;
         break;
      case 16:
         if(data_name!=NULL) sprintf(data_name,"Half-mass age [Gyrs]");
         if(data_type!=NULL) *data_type=SID_DOUBLE;
         if(data_d!=NULL)    *data_d   =fetch_treenode_delta_t_form(trees,current_halo)/S_PER_GYR;
         break;
      default:
         SID_trap_error("Invalid counter given in write_treenode_list_properties_set_ith()",ERROR_LOGIC);
   }
}

