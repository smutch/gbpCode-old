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

int fetch_treenode_n_particles_peak(tree_info *trees_in,tree_node_info *halo_in){
   if(halo_in!=NULL){
      // Use reference trees if present
      tree_info      *trees=fetch_trees_reference   (trees_in);
      tree_node_info *halo =fetch_treenode_reference(trees_in,halo_in);
      // Return the peak particle count
      if(halo!=NULL)
         return(halo->n_particles_peak);
   }
   return(-1);
}

