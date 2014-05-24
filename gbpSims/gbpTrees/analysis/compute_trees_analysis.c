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

void compute_trees_analysis(tree_info *trees){

  // Make sure the output directory exists
  mkdir(trees->filename_root_analysis,02755);

  // Compute merger tree statistics
  char filename_out_root[MAX_FILENAME_LENGTH];
  sprintf(filename_out_root,"%s/treenode",trees->filename_root_analysis);
  write_treenode_all_hist               (trees,filename_out_root);
  compute_trees_strayed_halo_analysis   (trees,filename_out_root);
  compute_trees_dropped_halo_analysis   (trees,filename_out_root);
  compute_trees_emerged_halo_analysis   (trees,filename_out_root);
  compute_trees_fragmented_halo_analysis(trees,filename_out_root);
  compute_trees_merger_analysis         (trees,filename_out_root);

}

