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

double fetch_treenode_delta_t(tree_info *trees,tree_node_info *halo_1,tree_node_info *halo_2){
   if(halo_1!=NULL && halo_2!=NULL){
      double t_1=trees->t_list[halo_1->snap_tree];
      double t_2=trees->t_list[halo_2->snap_tree];
      return(t_1-t_2);
   }
   return(0.);
}

