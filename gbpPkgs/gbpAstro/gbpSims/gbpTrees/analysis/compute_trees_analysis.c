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

  // Loop twice; once for groups and once for subgroups
  for(int i_type=0;i_type<2;i_type++){

     // Precompute marker statistics
     if(i_type==0)
        precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
     else
        precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);

     // Compute merger tree statistics
     write_treenode_all_hist                (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_analysis_strayed_halos   (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_analysis_dropped_halos   (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_analysis_emerged_halos   (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_analysis_fragmented_halos(trees,filename_out_root,i_type,logM_min,dlogM,n_logM);
     compute_trees_analysis_mergers         (trees,filename_out_root,i_type,logM_min,dlogM,n_logM);

     // Free precomputed markers
     if(i_type==0)
        free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_GROUPS);
     else
        free_precompute_treenode_markers(trees,PRECOMPUTE_TREENODE_MARKER_SUBGROUPS);
  }
}

