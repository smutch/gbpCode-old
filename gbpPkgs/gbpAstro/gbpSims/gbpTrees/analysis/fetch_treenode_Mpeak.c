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

double fetch_treenode_Mpeak(tree_info *trees_in,tree_node_info *halo_in){
   if(halo_in!=NULL){
      // Use reference trees if present
      tree_info      *trees=fetch_trees_reference   (trees_in);
      tree_node_info *halo =fetch_treenode_reference(trees_in,halo_in);
      // Fetch the peak mass marker
      tree_markers_info *markers  =fetch_treenode_precomputed_markers(trees,halo);
      if(markers!=NULL)
         return(markers->M_peak);
   }
   return(0.);
}

