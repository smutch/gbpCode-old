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

void compute_trees_analysis(tree_info *trees,double logM_min,double dlogM,int n_logM){

  // Make sure the output directory exists
  char filename_out_root[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root,"%s/treenode",trees->filename_root_analysis);
  mkdir(trees->filename_root_analysis,02755);

  for(int i_type=0;i_type<2;i_type++){

     // Precompute marker statistics
     if(i_type==0)
        precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
     else
        precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);

     // Compute merger tree statistics
     write_treenode_all_hist               (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_strayed_halo_analysis   (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_dropped_halo_analysis   (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_emerged_halo_analysis   (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_fragmented_halo_analysis(trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_merger_analysis         (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);

     if(i_type==0)
        free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
     else
        free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
  }
}

