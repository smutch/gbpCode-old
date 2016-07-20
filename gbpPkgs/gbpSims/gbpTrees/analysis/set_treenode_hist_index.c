#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

int set_treenode_hist_index(tree_info *trees,treenode_hist_info *hist,tree_node_info *current_halo,int i_axis){
   int     i_d=-1;
   int     i_prop;
   int    *args_i;
   double *args_d;
   int     flag_log;
   switch(i_axis){
      case 0:
        i_prop  =hist->x_prop;
        args_i  =hist->x_args_i;
        args_d  =hist->x_args_d;
        flag_log=hist->flag_log_x;
        break;
      case 1:
        i_prop  =hist->y_prop;
        args_i  =hist->y_args_i;
        args_d  =hist->y_args_d;
        flag_log=hist->flag_log_y;
        break;
   }
   // DEFINE THE BEHAVIOUR OF EACH PROPERTY HERE (don't forget to specify
   //   the elements of the treenode_hist_props structure, define the
   //   axis dimensions in init_treenode_hist() and describe the file
   //   writing in write_treenode_hist() as well)
   switch(i_prop){
      case 0:{ // z
        int di   =args_i[0];
        int i_min=trees->n_snaps%di;
        int i    =current_halo->snap_tree;
        i_d      =(i-i_min)/di;
        break;
      }
      case 1:{ // M_vir
        double d_min=args_d[0];
        double dd   =args_d[1];
        double d    =fetch_treenode_M_vir(trees,current_halo);
        if(flag_log)
           d=take_log10(d);
        i_d=(int)((d-d_min)/dd);
        break;
      }
      case 2:{ // M_peak
        double d_min=args_d[0];
        double dd   =args_d[1];
        double d    =fetch_treenode_M_peak(trees,current_halo);
        if(flag_log)
           d=take_log10(d);
        i_d=(int)((d-d_min)/dd);
        break;
      }
      case 3:{ // N
        double d_min=args_d[0];
        double dd   =args_d[1];
        double d    =(double)fetch_treenode_n_particles(trees,current_halo);
        if(flag_log)
           d=take_log10(d);
        i_d=(int)((d-d_min)/dd);
        break;
      }
      case 4:{ // N_peak
        double d_min=args_d[0];
        double dd   =args_d[1];
        double d    =(double)fetch_treenode_n_particles_peak(trees,current_halo);
        if(flag_log)
           d=take_log10(d);
        i_d=(int)((d-d_min)/dd);
        break;
      }
      case 5:{ // M_descendant
        double d_min=args_d[0];
        double dd   =args_d[1];
        double d    =fetch_treenode_M_vir(trees,current_halo->descendant);
        if(flag_log)
           d=take_log10(d);
        i_d=(int)((d-d_min)/dd);
        break;
      }
      case 6:{ // zeta=merger_ratio
        double d_min=args_d[0];
        double dd   =args_d[1];
        double d    =(double)fetch_treenode_zeta(trees,current_halo);
        if(flag_log)
           d=take_log10(d);
        i_d=(int)((d-d_min)/dd);
        break;
      }
      default:
        SID_trap_error("Invalid property passed to add_to_treenode_hist().",ERROR_LOGIC);
        break;
   }
   return(i_d);
}

