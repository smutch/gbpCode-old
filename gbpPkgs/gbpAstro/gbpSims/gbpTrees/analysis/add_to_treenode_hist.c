#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gbpLib.h>
#include <gbpMath.h>
#include <gbpHalos.h>
#include <gbpTrees_build.h>
#include <gbpTrees_analysis.h>
#include <assert.h>

void add_to_treenode_hist(tree_info *trees,treenode_hist_info *hist,tree_node_info *current_halo){
   int i_x=-1;
   int i_y=-1;
   int n_x=-1;
   int n_y=-1;
   for(int i_axis=0;i_axis<2;i_axis++){
      int    *i_d;
      int     i_prop;
      int    *args_i;
      double *args_d;
      int     flag_log;
      switch(i_axis){
         case 0:
           i_d     =&i_x;
           i_prop  =hist->x_prop;
           args_i  =hist->x_args_i;
           args_d  =hist->x_args_d;
           n_x     =hist->n_x;
           flag_log=hist->flag_log_x;
           break;
         case 1:
           i_d     =&i_y;
           i_prop  =hist->y_prop;
           args_i  =hist->y_args_i;
           args_d  =hist->y_args_d;
           n_y     =hist->n_y;
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
           (*i_d)   =(i-i_min)/di;
           break;
         }
         case 1:{ // M_vir
           double d_min=args_d[0];
           double dd   =args_d[1];
           double d    =fetch_treenode_Mvir(trees,current_halo);
           if(flag_log)
              d=take_log10(d);
           (*i_d)      =(int)((d-d_min)/dd);
           break;
         }
         case 2:{ // M_peak
           double d_min=args_d[0];
           double dd   =args_d[1];
           double d    =fetch_treenode_Mpeak(trees,current_halo);
           if(flag_log)
              d=take_log10(d);
           (*i_d)      =(int)((d-d_min)/dd);
           break;
         }
         case 3:{ // N
           double d_min=args_d[0];
           double dd   =args_d[1];
           double d    =(double)fetch_treenode_n_particles(trees,current_halo);
           if(flag_log)
              d=take_log10(d);
           (*i_d)      =(int)((d-d_min)/dd);
           break;
         }
         default:
           SID_trap_error("Invalid property passed to add_to_treenode_hist().",ERROR_LOGIC);
           break;
      }
   }

   // Populate the array 
   if(i_x>=0 && i_x<n_x && i_y>=0 && i_y<n_y)
      hist->array[i_y*n_x+i_x]++;
}

