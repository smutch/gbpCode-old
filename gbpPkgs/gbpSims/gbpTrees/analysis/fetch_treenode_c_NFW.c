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

double fetch_treenode_c_NFW(tree_info *trees,tree_node_info *halo){
   double c_NFW_halo=1.;
   static float **c_NFW_groups   =NULL;
   static float **c_NFW_subgroups=NULL;
   if(halo!=NULL){
      float **c_NFW;
      if(halo->parent==NULL){
         if(c_NFW_groups==NULL)
            c_NFW_groups=(float **)ADaPS_fetch(trees->data,"c_NFW_groups");
         c_NFW=c_NFW_groups;
      }
      else{
         if(c_NFW_subgroups==NULL)
            c_NFW_subgroups=(float **)ADaPS_fetch(trees->data,"c_NFW_subgroups");
         c_NFW=c_NFW_subgroups;
      }
      float *c_NFW_snap=c_NFW[fetch_treenode_snap_tree(trees,halo)];
      if(c_NFW_snap!=NULL)
         c_NFW_halo=c_NFW_snap[fetch_treenode_file_index(trees,halo)];
      else
         SID_trap_error("NFW concentrations have not been read for snapshot %d.",trees->snap_list[fetch_treenode_snap_tree(trees,halo)]);
   }
   return((double)c_NFW_halo);
}

