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

tree_markers_info *fetch_treenode_precomputed_markers(tree_info *trees,tree_node_info *halo){
   tree_markers_info *markers=NULL;
   if(halo!=NULL){
      tree_markers_info **markers_precomputed=NULL;
      if(halo->parent!=NULL)
         markers_precomputed=trees->subgroup_markers;
      else
         markers_precomputed=trees->group_markers;
      if(markers_precomputed==NULL)
         SID_trap_error("Precomputed markers have not been made available.",ERROR_LOGIC);
      markers=&(markers_precomputed[halo->snap_tree][halo->neighbour_index]);
   }
   if(markers==NULL)
      SID_trap_error("Could not set the markers for a halo in fetch_treenode_precomputed_markers(). (is_halo_NULL=%d)",ERROR_LOGIC,halo==NULL);
   return(markers);
}

