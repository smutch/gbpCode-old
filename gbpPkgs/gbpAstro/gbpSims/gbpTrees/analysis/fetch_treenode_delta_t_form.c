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

double fetch_treenode_delta_t_form(tree_info *trees,tree_node_info *halo){
   if(halo!=NULL){
      tree_markers_info *markers=fetch_treenode_precomputed_markers(trees,halo);
      if(markers!=NULL){
         tree_node_info *form=markers->half_peak_mass;
         if(form!=NULL){
            double t_1=trees->t_list[halo->snap_tree];
            double t_2=trees->t_list[form->snap_tree];
            return(t_1-t_2);
         }
      }
   }
   return(0.);
}

