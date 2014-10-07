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

int find_treenode_markers(tree_info *trees,tree_node_info *halo,tree_markers_info **markers_all,tree_markers_info *markers){
   if(markers_all==NULL){
      find_treenode_branch_root    (trees,halo,&(markers->branch_root));
      find_treenode_branch_leaf    (trees,halo,&(markers->branch_leaf));
      find_treenode_main_progenitor(trees,halo,&(markers->main_progenitor));
      find_treenode_accretion      (trees,halo,&(markers->first_became_satellite),
                                               &(markers->joined_current_parent));
      find_treenode_formation      (trees,halo,0.5,
                                               &(markers->peak_mass),
                                               &(markers->half_peak_mass));
      find_treenode_last_merger    (trees,halo,0.33,
                                               &(markers->merger_33pc_remnant),
                                               &(markers->merger_33pc_host),
                                               &(markers->merger_33pc_merger));
      find_treenode_last_merger    (trees,halo,0.1,
                                               &(markers->merger_10pc_remnant),
                                               &(markers->merger_10pc_host),
                                               &(markers->merger_10pc_merger));
      if(markers->branch_root!=NULL)
         markers->descendant=markers->branch_root->descendant;
      else
         markers->descendant=NULL;
      markers->flag_halo_is_main_progenitor=check_treenode_if_main_progenitor(halo);
   }
   else
      memcpy(markers,&(markers_all[halo->snap_tree][halo->neighbour_index]),sizeof(tree_markers_info));
}

